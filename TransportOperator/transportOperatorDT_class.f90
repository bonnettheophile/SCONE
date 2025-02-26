!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng
  use rng_class,                  only : rng

  ! Superclass
!  use transportOperator_inter,    only : transportOperator
  use transportOperator_inter,    only : transportOperator, init_super => init     ! Added init routine to enable case runs directly from input files
!  use transportOperator_inter,    only : transportOperator
  use transportOperator_inter,    only : transportOperator, init_super => init     ! Added init routine to enable case runs directly from input files

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,      only : nucReg_get => get
  use nuclearDatabase_inter,       only : nuclearDatabase
  use materialMenu_mod,            only : mm_matIdx => matIdx
  
  use nuclearDataReg_mod,      only : nucReg_get => get
  use nuclearDatabase_inter,       only : nuclearDatabase
  use materialMenu_mod,            only : mm_matIdx => matIdx
  
  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator):: transportOperatorDT
!   Data definitions for virtual density module  !
    real(defReal)                    :: product_factor 
    real(defReal), dimension(3)      :: vector_factor, vector_factor_cur
    logical(defBool)                 :: virtual_density, cross_over = .false.
    character(nameLen)               :: deform_type, direction_type, scale_type
    character(nameLen), allocatable  :: pert_mat(:)
    integer(shortInt), allocatable   :: pert_mat_id(:)
    integer(shortInt)                :: nb_pert_mat
    logical(defBool)                 :: swelled
!   Data definitions for virtual density module  !
  contains
    procedure:: transit => deltaTracking
!   Virtual density addition  !
    procedure:: init
!   Virtual density addition  !

  end type transportOperatorDT

contains

  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout):: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter:: Here = 'deltaTracking (transportOIperatorDT_class.f90)'
    character(100):: scale
!   Data definitions for virtual density module  !
    real(defReal), dimension(3)                :: cosines, virtual_cosines, real_vector, virtual_vector
    real(defReal)                             :: virtual_dist, flight_stretch_factor

!   Data definitions for virtual density end here!

    ! Get majorat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE/self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)
    scale = trim(self % scale_type)

    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p % pRNG % get() ) * majorant_inv*self % product_factor
      if (self % virtual_density) then
        cosines(:) = p % dirGlobal()
        real_vector = distance*cosines

        if (self % deform_type == 'swelling') then
          virtual_vector(1) = real_vector(1) * self % vector_factor(2) * self % vector_factor(3)
          virtual_vector(2) = real_vector(2) * self % vector_factor(1) * self % vector_factor(3)
          virtual_vector(3) = real_vector(3) * self % vector_factor(1) * self % vector_factor(2)
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist/distance
          virtual_cosines(1) = cosines(1) * self % vector_factor(2) * self % vector_factor(3) / flight_stretch_factor
          virtual_cosines(2) = cosines(2) * self % vector_factor(1) * self % vector_factor(3) / flight_stretch_factor
          virtual_cosines(3) = cosines(3) * self % vector_factor(1) * self % vector_factor(2) / flight_stretch_factor
        elseif (self % deform_type == 'expansion') then
          virtual_vector = real_vector/self % vector_factor
          virtual_dist = sqrt(sum(virtual_vector**2))
          flight_stretch_factor = virtual_dist/distance
          virtual_cosines = cosines / (self % vector_factor*flight_stretch_factor)
        else
          print *,'Error in recognizing type of geometric deformation  ! Please check input!'
        end if
      
        call p % point(virtual_cosines)
        distance = virtual_dist
      end if
    ! Move partice in the geometry

      call self % geom % teleport(p % coords, distance)         

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      ! Check for void
      if (p % matIdx() == VOID_MAT) then
        call tally % reportInColl(p, .true.)
        cycle DTLoop
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

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
    class(transportOperatorDT), intent(inout):: self
    class(dictionary), intent(in)             :: dict
    real(defReal), allocatable                :: vec(:)
    character(100), parameter                 :: Here = 'init (transportOperatorDT_class.f90)'
         !Virtual Density Data call  begins !

    call init_super(self, dict)
    call dict % getorDefault(self % virtual_density, 'virtual_density', .false.)
    call dict % getorDefault(self % swelled, 'swelled', .false.)
    self % product_factor = ONE
    if (self % swelled) then
        call dict % get(vec, "factor_1")
        self % vector_factor = vec

        self % product_factor = product(self % vector_factor)
    end if

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
