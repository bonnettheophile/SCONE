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

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
!   Data definitions for virtual density module  !
    real(defReal)                    :: product_factor 
    real(defReal), dimension(3)      :: vector_factor, vector_factor_cur
    logical(defBool)                 :: virtual_density, cross_over = .false.
    character(nameLen)               :: deform_type, direction_type, scale_type
    character(nameLen), allocatable  :: pert_mat(:)
    integer(shortInt), allocatable   :: pert_mat_id(:)
    integer(shortInt)                :: nb_pert_mat
!   Data definitions for virtual density module  !
  contains
    procedure :: transit => deltaTracking
!   Virtual density addition  !
    procedure :: init
!   Virtual density addition  !

  end type transportOperatorDT

contains

  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance, speed, time
    character(100), parameter :: Here = 'deltaTracking (transportOIperatorDT_class.f90)'
    character(100) :: scale
!   Data definitions for virtual density module  !
    real(defReal),dimension(3)                :: cosines,virtual_cosines, real_vector, virtual_vector, preCosines
    real(defReal)                             :: virtual_dist, flight_stretch_factor

!   Data definitions for virtual density end here!

    ! Get majorat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)
    preCosines = p % dirGlobal()

    scale = trim(self % scale_type)

    ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p % pRNG % get() ) * majorant_inv
      if (self % virtual_density) then
        cosines(:) = preCosines
        real_vector = distance * cosines

        if (self % deform_type == 'swelling') then
          virtual_vector(1) = real_vector(1) * self % vector_factor(2) * self % vector_factor(3)
          virtual_vector(2) = real_vector(2) * self % vector_factor(1) * self % vector_factor(3)
          virtual_vector(3) = real_vector(3) * self % vector_factor(1) * self % vector_factor(2)
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist / distance
          virtual_cosines(1) = cosines(1) * self % vector_factor(2) * self % vector_factor(3) / flight_stretch_factor
          virtual_cosines(2) = cosines(2) * self % vector_factor(1) * self % vector_factor(3) / flight_stretch_factor
          virtual_cosines(3) = cosines(3) * self % vector_factor(1) * self % vector_factor(2) / flight_stretch_factor
        elseif (self % deform_type == 'expansion') then
          virtual_vector = real_vector / self % vector_factor
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist/distance
          virtual_cosines = cosines / (self % vector_factor*flight_stretch_factor)
        else
          print *,'Error in recognizing type of geometric deformation! Please check input!'
        end if
      
        call p % point(virtual_cosines)
        distance = virtual_dist
      end if
    ! Move partice in the geometry

      call self % geom % teleport(p % coords, distance)   
      call p % point(preCosines)  

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

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

  subroutine init(self, dict)
    class(transportOperatorDT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    real(defReal), allocatable                :: vec(:)
    character(100), parameter                 :: Here = 'init (transportOperatorDT_class.f90)'
         !Virtual Density Data call  begins !

    call init_super(self, dict)
    call dict % getorDefault(self % virtual_density, 'virtual_density', .false.)
      if (self % virtual_density) then
        call dict % getorDefault(self % deform_type, 'deform_type_1','swelling')
        call dict % getorDefault(self % direction_type, 'direction_type','isotropic')
        call dict % getorDefault(self % scale_type, 'scale','uniform')
        call dict % get(vec, "factor_1")
        self % vector_factor = vec

        self % product_factor = product(self % vector_factor)

        if (trim(self % scale_type) == 'non_uniform') then
          call fatalError(Here, "Delta tracking cannot be used with non-uniform virtual densities")
        end if
      end if
!    Virtual Density Data call  ends !
  end subroutine init
  
end module transportOperatorDT_class
