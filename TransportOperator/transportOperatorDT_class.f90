!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use errors_mod,               only : fatalError
  use genericProcedures,        only : numToChar
  use particle_class,           only : particle
  use particleDungeon_class,    only : particleDungeon
  use dictionary_class,         only : dictionary

  ! Superclass
  use transportOperator_inter,  only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,           only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,         only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,       only : ndReg_get => get
  use nuclearDatabase_inter,    only : nuclearDatabase
  use materialMenu_mod,            only : mm_matIdx => matIdx


  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
    real(defReal)                    :: product_factor 
    real(defReal), dimension(:,:), allocatable      :: vector_factor
    logical(defBool)                 :: virtual_density, cross_over = .false.
    character(nameLen)               :: direction_type, scale_type
    character(nameLen), allocatable  :: pert_mat(:), deform_type(:)
    integer(shortInt), allocatable   :: pert_mat_id(:)
    integer(shortInt)                :: nb_pert_mat
  contains
    procedure :: transit => deltaTracking
    ! Override procedure
    procedure :: init

  end type transportOperatorDT

contains

  !!
  !! Performs delta tracking until a real collision point is found
  !!
  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance, speed, time
    real(defReal)                             :: cosines(3), real_vector(3), virtual_vector(3)
    real(defReal)                             :: virtual_cosines(3), virtual_dist, flight_stretch_factor
    character(100), parameter :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! Get majorant XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      call p % pointSaved(p % dirGlobal())
      distance = -log( p% pRNG % get() ) * majorant_inv
        
      speed = p % getSpeed()
      time = distance / speed + p % time

      ! Set a max flight distance due to hitting the time-boundary
      if (p % timeMax > ZERO .and. time > p % timeMax) then
        distance = speed * (p % timeMax - p % time)
        p % fate = AGED_FATE
      end if

      if (self % virtual_density .and. trim(self % scale_type) == 'uniform') then
        cosines(:) = p % coords % lvl(1) % savedDir 
        real_vector = distance*cosines

        if (self % deform_type(1) == 'swelling') then
          virtual_vector(1) = real_vector(1) * p % f(2) * p % f(3)
          virtual_vector(2) = real_vector(2) * p % f(1) * p % f(3)
          virtual_vector(3) = real_vector(3) * p % f(1) * p % f(2)
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist/distance
          virtual_cosines(1) = cosines(1) * p % f(2) * &
              p % f(3) / flight_stretch_factor
          virtual_cosines(2) = cosines(2) * p % f(1) * &
              p % f(3) / flight_stretch_factor
          virtual_cosines(3) = cosines(3) * p % f(1) * &
              p % f(2) / flight_stretch_factor
        elseif (self % deform_type(1) == 'expansion') then
          virtual_vector = real_vector/p % f
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist/distance
          virtual_cosines = cosines / (p % f*flight_stretch_factor)
        else
          call fatalError(Here, 'Unrecognised geometric deformation')
        end if
      
        call p % point(virtual_cosines)

        distance = virtual_dist
      end if

      ! Move particle in the geometry and time
      ! Note: teleport routine has been modified to use savedDir instead of current dir
      call self % geom % teleport(p % coords, distance)
      p % time = p % time + distance / speed
      
      select case(p % matIdx())

        ! If particle has leaked exit
        case(OUTSIDE_FILL)
          p % fate = LEAK_FATE
          p % isDead = .true.
          exit DTLoop

        ! Check for void
        case(VOID_MAT)
          if (p % fate == AGED_FATE) exit DTLoop
          call tally % reportInColl(p, .true.)
          cycle DTLoop

        ! Give error if the particle somehow ended in an undefined material
        case(UNDEF_MAT)
          print*, 'Particle location: ', p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")

        ! Give error if the particle is in a region with overlapping cells
        case(OVERLAP_MAT)
          print*, 'Particle location: ', p % rGlobal()
          call fatalError(Here, "Particle is in overlapping cells")

        case default
          ! All is well        

      end select
      
      ! If particle has aged, exit
      if (p % fate == AGED_FATE) then
        exit DTLoop
      end if

      ! Get local conditions of temperature and density
      call self % localConditions(p)
      
      ! Obtain the local cross-section
      sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (p % pRNG % get() < sigmaT*majorant_inv) then
        exit DTLoop
      else
        call tally % reportInColl(p, .true.)
      end if

    end do DTLoop

    call tally % reportTrans(p)

  end subroutine deltaTracking

  !!
  !! Initialise DT transport operator
  !!
  !! See transportOperator_inter for more details
  !!
  subroutine init(self, dict)
    class(transportOperatorDT), intent(inout):: self
    class(dictionary), intent(in)             :: dict
    character(nameLen)                        :: input
    real(defReal), allocatable, dimension(:)  :: vec                           
    integer(shortInt)                         :: index

    ! Initialise superclass
    call init_super(self, dict)

    ! Initialise virtual density
    call dict % getorDefault(self % virtual_density, 'virtual_density', .false.)
    if (self % virtual_density) then
      call dict % getorDefault(self % direction_type, 'direction_type','isotropic')
      call dict % getorDefault(self % scale_type, 'scale','uniform')
      call dict % getorDefault(self % nb_pert_mat, 'nb_pert_mat', 1)

      if (trim(self % scale_type) == 'non_uniform') then
        allocate(self % deform_type(self % nb_pert_mat))
        allocate(self % pert_mat(self % nb_pert_mat))
        allocate(self % pert_mat_id(self % nb_pert_mat))
        allocate(self % vector_factor(3, self % nb_pert_mat))
        do index = 1, self % nb_pert_mat
          input = 'factor_'
          write(input, '(I0)') index
          input = trim('factor_')//trim(input)
          call dict % get(vec, trim(input))
          self % vector_factor(:,index) = vec

          input = 'pert_mat_'
          write(input, '(I0)') index
          input = trim('pert_mat_')//trim(input)
          call dict % get(self % pert_mat(index), trim(input))
          self % pert_mat_id(index) = mm_matIdx(self % pert_mat(index))

          input = 'deform_type_'
          write(input, '(I0)') index
          input = trim('deform_type_')//trim(input)
          call dict % get(self % deform_type(index), trim(input))
        end do
      else
        allocate(self % deform_type(1))
        allocate(self % vector_factor(3, 1))
        call dict % get(self % deform_type(1), "deform_type_1")
        call dict % get(vec, "factor_1")
        self % vector_factor(:,1) = vec
      end if
    end if

  end subroutine init


end module transportOperatorDT_class
