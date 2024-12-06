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
    character(nameLen)                        :: test_name = 'clad'
    real(defReal)                             :: enquiry_XS_test
    real(defReal)                             :: enquiry_XS_max
    integer(shortInt)                         :: mat_Idx_unpert = 0, mat_Idx_pert = 0
    logical(defBool)                          :: weight_pert, denizen_pert, pert_region

!   Data definitions for virtual density end here!

    ! Get majorat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p % pRNG % get() ) * majorant_inv
       if (self % virtual_density) then
          if (trim(self % scale_type) == 'non_uniform') then
               if (any(self % pert_Mat_id == p % matIdx())) then
                self % vector_factor = self % vector_factor_cur
                self % deform_type = 'swelling'
               else
                 self % vector_factor = ONE
               end if
           end if

          if (((trim(self % scale_type) == 'non_uniform') .and. (any(self % pert_Mat_id == p % matIdx()))) &
              .or. (trim(self % scale_type) == 'uniform')) then    
                cosines(:) = p % dirGlobal()
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

               denizen_pert = .true.    ! Indicates that this neutron has been found inside perturbed region
           end if
          if (trim(self % scale_type) == 'non_uniform') then
            if (any(self % pert_Mat_id == p % matIdx())) then
               pert_region = .true.         ! indicates neutron is in perturbed region
               mat_Idx_pert = p % matIdx()
               mat_Idx_unpert = 0
             else
               pert_region = .false.         ! indicates prisoner is in unperturbed region
               mat_Idx_unpert = p % matIdx()
               mat_Idx_pert = 0
            end if
         end if
       end if
      ! Move partice in the geometry

      call self % geom % teleport(p % coords, distance)

      if ((self % virtual_density) .and. (trim(self % scale_type) == 'non_uniform')) then

         if ((pert_region) .and. (all(self % pert_Mat_id /= p % matIdx()))) then !&.and. (p % matIdx() /= OUTSIDE_FILL)) then
              pert_region = .false.

              if ((abs(p % w - ONE) < 1e-5) .or. (abs(p % w - TWO) < 1e-5)) then ! prisoner is crossing over from its original perturbed cell to an unperturbed region
                self % cross_over = .true.         ! marks the prisoner to have crossed over to a new region
                p % w = p % w / (1.0 + ((self % vector_factor(3) - ONE) / ONE))   ! Weight is perturbed on cross-over to a new region  
              elseif ((abs(p % w - (ONE / (ONE + ((self % vector_factor(3) - ONE) / ONE)))) < 1e-5) &
               .or. (abs(p % w - (TWO / (ONE + ((self % vector_factor(3) - ONE) / ONE)))) < 1e-5)) then ! prisoner has now gone back to his original region!
                self % cross_over = .true.         ! markes the prisoner to have crossed over to a new region !
                p % w = p % w * (ONE + ((self % vector_factor(3) - ONE) / ONE))  ! Weight is perturbed on cross-over to a new region
              end if

          elseif  ((.not. pert_region) .and. (any(self % pert_Mat_id == p % matIdx()))) then
              pert_region = .true.
              if ((abs(p % w - ONE) < 1e-5) .or. (abs(p % w - TWO) < 1e-5)) then ! prisoner is crossing over from his original region to a new region !
                self % cross_over = .true.         ! markes the prisoner to have crossed over to a new region !
                p % w = p % w / (ONE + ((self % vector_factor(3) - ONE) / ONE))   ! Weight is perturbed on cross-over to a new region
              elseif ((abs(p % w - (ONE / (ONE + ((self % vector_factor(3) - ONE) / ONE)))) < 1e-5) &
               .or. (abs(p % w - (TWO / (ONE + ((self % vector_factor(3) - ONE) / ONE)))) < 1e-5)) then ! prisoner has now gone back to his original region!
                self % cross_over = .false.         ! markes the prisoner is now back in his original region !
                p % w = p % w * (ONE + ((self % vector_factor(3) - ONE) / ONE))    ! Weight perturbation is removed since the prisoner no longer contributes to the current
             end if
                  
          end if
            if ((p % w < 0.88) .and. ((self % vector_factor(3) - ONE) > 1e-5)) then
                call fatalError(Here, "Warning! Incorrect current perturbation!")
                print *,'mat_Idx_unpert,mat_Idx_pert, p % matIdx() = ',mat_Idx_unpert, mat_Idx_pert, p % matIdx()
                print *,'p % w = ',p % w
                print *,'pert_region_cond, self % cross_over_cond = ', pert_region, self % cross_over
            end if                
       end if

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
            !if ((abs(p % w - ONE / (ONE + ((self % vector_factor(3) - ONE) / TWO))) < 1e-5) .or. &
             !(abs(p % w - TWO / (ONE + ((self % vector_factor(3) - ONE) / TWO))) < 1e-5)) then ! prisoner has now gone back to his original region!
                !p % w = p % w * (ONE + ((self % vector_factor(3) - ONE) / TWO))
                !print *,'p % w = ',p % w
            !end if
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
