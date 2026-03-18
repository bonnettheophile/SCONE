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
  use geometryReg_mod,                only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use deformationField_class,         only : deformationField, deformationField_TptrCast


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
    real(defReal), dimension(:), allocatable        :: vectorFactor
    logical(defBool)                                :: virtualDensity
    logical(defBool)                                :: isGPC
    character(nameLen)                              :: scaleType
    character(nameLen), allocatable                 :: deformType(:)

    type(deformationField), pointer :: deformationField
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
    real(defReal)                             :: cosines(3), realVector(3), virtualVector(3)
    real(defReal)                             :: virtualCosines(3), virtualDist, flightStretch
    real(defReal)                             :: f(3)
    character(100), parameter :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! If gpc calculation, use particle-wise deformation
    if (self % isGPC .and. trim(self % scaleType) == 'uniform') then
      f = p % f
    elseif (trim(self % scaleType) == 'uniform') then
      f = self % vectorFactor
    end if

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

      if (self % virtualDensity .and. trim(self % scaleType) == 'uniform') then
        cosines(:) = p % coords % lvl(1) % savedDir 
        realVector = distance*cosines

        if (self % deformType(1) == 'swelling') then
          virtualVector(1) = realVector(1) * f(2) * f(3)
          virtualVector(2) = realVector(2) * f(1) * f(3)
          virtualVector(3) = realVector(3) * f(1) * f(2)
          virtualDist = sqrt(sum(virtualVector**2))
          flightStretch = virtualDist/distance
          virtualCosines(1) = cosines(1) * f(2) * &
              f(3) / flightStretch
          virtualCosines(2) = cosines(2) * f(1) * &
              f(3) / flightStretch
          virtualCosines(3) = cosines(3) * f(1) * &
              f(2) / flightStretch
        elseif (self % deformType(1) == 'expansion') then
          virtualVector = realVector / f
          virtualDist = sqrt(sum(virtualVector**2))
          flightStretch = virtualDist / distance
          virtualCosines = cosines / (f * flightStretch)
        else
          call fatalError(Here, 'Unrecognised geometric deformation')
        end if
      
        call p % point(virtualCosines)

        distance = virtualDist

      elseif (self % virtualDensity .and. trim(self % scaleType) == 'non_uniform') then

        ! Transform to real space
        if (self % isGPC) then
          realVector = self % deformationField % forward(p % coords, p % X)
        else
          realVector = self % deformationField % forward(p % coords)
        end if

        call p % coords % assignPosition(realVector)
      end if

      ! Move particle in the geometry and time
      ! Note: teleport routine has been modified to use savedDir instead of current dir
      ! Note: boundary conditions should not be perturbed by local deformation
      call self % geom % teleport(p % coords, distance)

      if (self % virtualDensity .and. trim(self % scaleType) == 'non_uniform') then

        ! Transform back to material space for xs lookup
        if (self % isGPC) then
          realVector = self % deformationField % backward(p % coords, p % X)
        else
          realVector = self % deformationField % backward(p % coords)
        end if

        call p % coords % assignPosition(realVector)
        call self % geom % placeCoord(p % coords)
      end if

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
    real(defReal), allocatable, dimension(:)  :: vec                           
    integer(shortInt)                         :: index

    ! Initialise superclass
    call init_super(self, dict)

    ! Initialise virtual density
    call dict % getorDefault(self % virtualDensity, 'virtualDensity', .false.)
    if (self % virtualDensity) then
      call dict % getorDefault(self % scaleType, 'scale','uniform')
      call dict % getorDefault(self % isGPC, 'gpc', .false.)

      if (trim(self % scaleType) == 'non_uniform') then
        index = gr_fieldIdx(nameDeformation)
        self % deformationField => deformationField_TptrCast(gr_fieldPtr(index))
      else
        allocate(self % deformType(1))
        allocate(self % vectorFactor(3))
        call dict % get(self % deformType(1), "deformType")
        call dict % get(vec, "factor")
        self % vectorFactor = vec
      end if
    end if

  end subroutine init


end module transportOperatorDT_class
