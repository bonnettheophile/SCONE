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

  ! Superclass
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
  
  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
!   Data definitions for virtual density module  !
    real(defReal)                  :: product_factor 
    real(defReal), dimension(3)    :: vector_factor, vector_factor_cur
    logical(defBool)               :: virtual_density, cross_over = .false.
    character(nameLen)             :: deform_type, direction_type, scale_type
    character(nameLen)             :: perturb_mat
    integer(shortInt),dimension(6) :: pert_Mat_Id = 0
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
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOIperatorDT_class.f90)'
!   Data definitions for virtual density module  !
    real(defReal),dimension(3)                :: cosines,virtual_cosines, real_vector, virtual_vector, virtual_XS
    real(defReal)                             :: virtual_dist, flight_stretch_factor, virtual_XS_max_inv
    character(nameLen)                        :: test_name = 'fuel', test_name2 = 'clad'
    real(defReal),dimension(2)                             :: enquiry_XS_test
    real(defReal)                             :: enquiry_XS_max
    integer(shortInt)                         :: mat_Idx_unpert = 0, mat_Idx_pert = 0
    logical(defBool)                          :: weight_pert, denizen_pert, pert_region

!   Data definitions for virtual density end here!

    ! Get majorat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    if (self % virtual_density .and. .false.) then   ! Preliminary check to see if majorant needs to be checked against virtual cross-sections!

          if (((self % deform_type == 'swelling') .and. (self % product_factor < ONE )) .or. &
             (((self % deform_type == 'expansion') .and. (self % product_factor > ONE )))) then
            enquiry_XS_test(1) = self % xsData % getTrackMatXS(p, mm_matidx(test_name))
            enquiry_XS_test(2) = self % xsData % getTrackMatXS(p, mm_matidx(test_name2))
            enquiry_XS_max=maxval(enquiry_XS_test)
            !enquiry_XS_max = enquiry_XS_test ! because not an array in this case

              if (self % deform_type == 'swelling') then
               virtual_XS(1) = enquiry_XS_max / (self % vector_factor(2) * self % vector_factor(3))
               virtual_XS(2) = enquiry_XS_max / (self % vector_factor(1) * self % vector_factor(3))
               virtual_XS(3) = enquiry_XS_max / (self % vector_factor(1) * self % vector_factor(2))
              elseif (self % deform_type == 'expansion') then
               virtual_XS = enquiry_XS_max * self % vector_factor
              end if
              virtual_XS_max_inv= ONE / maxval(virtual_XS)      ! Inverse of the maximum virtual cross-section!
              if (virtual_XS_max_inv < majorant_inv) then    ! Condition to check if the majorant will be exceeded in virtual perturbation!!
                 majorant_inv = virtual_XS_max_inv
              end if
          end if
    end if

    DTLoop:do
      distance = -log( p % pRNG % get() ) * majorant_inv
      if (self % virtual_density) then
        if (trim(self % scale_type) == 'non_uniform') then
          if (any(self % pert_Mat_id == p % matIdx())) then
            p % isPerturbed = .true. ! Set particle to be perturbed
            if (p % startPerturbed == 0) then
              p % startPerturbed = 1 ! Set particle starting in perturbed region
              p % lastPerturbed = .true.
            end if
          else
            p % isPerturbed = .false.
            if (p % startPerturbed == 0) then 
              p % startPerturbed = 2 ! Set particle starting in unperturbed region
              p % lastPerturbed = .false.
            end if
          end if
        end if

          if (((trim(self % scale_type) == 'non_uniform') .and. (any(self % pert_Mat_id == p % matIdx()))) &
              .or. (trim(self % scale_type) == 'uniform')) then    
                cosines(:) = p % dirGlobal()
                real_vector = distance * cosines
                
               if (self % deform_type == 'swelling') then
                  virtual_vector(1) = real_vector(1) * p % f(2) * p % f(3)
                  virtual_vector(2) = real_vector(2) * p % f(1) * p % f(3)
                  virtual_vector(3) = real_vector(3) * p % f(1) * p % f(2)
                  virtual_dist = sqrt(sum(virtual_vector**2))
                  flight_stretch_factor = virtual_dist / distance
                  virtual_cosines(1) = cosines(1) * p % f(2) * p % f(3) / flight_stretch_factor
                  virtual_cosines(2) = cosines(2) * p % f(1) * p % f(3) / flight_stretch_factor
                  virtual_cosines(3) = cosines(3) * p % f(1) * p % f(2) / flight_stretch_factor
               elseif (self % deform_type == 'expansion') then
                  virtual_vector = real_vector / p % f
                  virtual_dist = sqrt(sum(virtual_vector**2))
                  flight_stretch_factor = virtual_dist/distance
                  virtual_cosines = cosines / (p % f * flight_stretch_factor)
               else
                  print *,'Error in recognizing type of geometric deformation! Please check input!'
               end if
              
               call p % point(virtual_cosines)
               distance = virtual_dist

        end if
      end if
    ! Move partice in the geometry

      call self % geom % teleport(p % coords, distance)

      if (trim(self % scale_type) == 'non_uniform') then
        if (any(self % pert_Mat_id == p % matIdx())) then
            p % isPerturbed = .true.
          else
            p % isPerturbed = .false.
        end if
      end if

      if ((self % virtual_density) .and. (trim(self % scale_type) == 'non_uniform')) then
        if (( p % lastPerturbed .eqv. .false.) .and. p % isPerturbed) then
          if (p % startPerturbed == 1) then
            !p % w = p % w / self % vector_factor(3)
            p % w = p % w / self % vector_factor(3)**2 !isotropic swelling
          else
            p % w = p % w * self % vector_factor(3)**2 !isotropic swelling
            !p % w = p % w * self % vector_factor(3)
          end if
        elseif (( p % lastPerturbed) .and. (p % isPerturbed .eqv. .false.)) then
          if (p % startPerturbed == 1) then
            p % w = p % w * self % vector_factor(3)**2 !isotropic swelling
            !p % w = p % w / self % vector_factor(3)
          else
            p % w = p % w / self % vector_factor(3)**2 !isotropic swelling
            !p % w = p % w * self % vector_factor(3)
          end if
        end if
        p % lastPerturbed = p % isPerturbed
      end if              

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
    class(transportOperatorDT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
         !Virtual Density Data call  begins !

    call init_super(self, dict)
    call dict % getorDefault(self % virtual_density, 'virtual_density', .false.)
      if (self % virtual_density) then
        call dict % getorDefault(self % deform_type, 'deform_type','swelling')
        call dict % getorDefault(self % direction_type, 'direction_type','isotropic')
        call dict % getorDefault(self % scale_type, 'scale','uniform')

        if (trim(self % direction_type) == 'anisotropic') then
          call dict % getorDefault(self % vector_factor(1), 'x_factor', ONE)
          call dict % getorDefault(self % vector_factor(2), 'y_factor', ONE)
          call dict % getorDefault(self % vector_factor(3), 'z_factor', ONE)
        else
          call dict % getorDefault(self % vector_factor(1), 'factor', ONE)
          call dict % getorDefault(self % vector_factor(2), 'factor', ONE)
          call dict % getorDefault(self % vector_factor(3), 'factor', ONE)
        end if

        self % vector_factor_cur = self % vector_factor
        self % product_factor = self % vector_factor(1) * self % vector_factor(2) * self % vector_factor(3)

        if (trim(self % scale_type) == 'non_uniform') then
              call dict % getorDefault(self % perturb_mat, 'pert_mat','uniform')
              self % pert_Mat_Id(1) = mm_matIdx(self % perturb_mat)
        end if
      end if
!    Virtual Density Data call  ends !
  end subroutine init
  
end module transportOperatorDT_class
