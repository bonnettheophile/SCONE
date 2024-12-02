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
    real(defReal)                  :: x_factor,y_factor,z_factor,deform_factor
    integer(shortInt)              :: Vir_Den_Swch,cross_over_cond = 0
    character(nameLen)             :: deform_type, direction_type, scale_type
    character(nameLen)             :: perturb_mat1,perturb_mat2, perturb_mat3,perturb_mat4,perturb_mat5
    character(nameLen)             :: perturb_mat6,perturb_mat7,perturb_mat8,perturb_mat9,perturb_mat10
    character(nameLen)             :: perturb_mat11,perturb_mat12
    real(defReal)                  :: x_factor_cur,y_factor_cur,z_factor_cur
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
    real(defReal),dimension(3)                :: cosines,virtual_cosines,dist_vector,Vir_Den_dist_vector,virtual_XS
    real(defReal)                             :: Vir_Den_dist,flight_stretch_factor,virtual_XS_max_inv
    character(nameLen)                        :: test_name1 = 'clad',test_name2 = 'can',test_name3 = 'lead',test_name4 = 'gap'
    character(nameLen)                        :: test_name5 = 'fuel',test_name6 = 'b4c'
    real(defReal),dimension(6)                :: enquiry_XS_test
    real(defReal)                             :: enquiry_XS_max
    integer(shortInt)                         :: pert_region_cond = 0, mat_Idx_unpert = 0, mat_Idx_pert = 0, denizen_pert
    integer(shortInt)                         :: weight_pert_cond

!   Data definitions for virtual density end here!

    ! Get majorat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    if (self % Vir_Den_Swch==1) then   ! Preliminary check to see if majorant needs to be checked against virtual cross-sections!
        denizen_pert = 0         ! Initialization of perturbation identifier! 0 indicates that at the start of the loop, this neutron is not perturbed!
        weight_pert_cond = 0     ! Initialization of perturbation identifier! 0 indicates that at the start of the loop, the weight of this neutron is not perturbed!
        pert_region_cond = 0     ! Initialization of perturbed region condition as unperturbed! This is prior to start of loop!
        !print *,'Top denizen_pert, weight_pert_cond, pert_region_cond = ', denizen_pert, weight_pert_cond, pert_region_cond   ! Confirmed that these parameters are correctly initialized!

          if (((self % deform_type == 'swelling') .and. (self % x_factor * self % y_factor * self % z_factor < 1.0 )) .or. &
             (((self % deform_type == 'expansion') .and. (self % x_factor * self % y_factor * self % z_factor > 1.0 )))) then
            enquiry_XS_test(1) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name1))
            enquiry_XS_test(2) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name2))
            enquiry_XS_test(3) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name3))
            enquiry_XS_test(4) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name4))
            enquiry_XS_test(5) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name5))
            enquiry_XS_test(6) = self % xsData % getTotalMatXS(p, mm_matIdx(test_name6))
            enquiry_XS_max=maxval(enquiry_XS_test)
            !print *,'mm_matIdx(test_name1) = ',mm_matIdx(test_name1)
            !print *,'mm_matIdx(test_name2) = ',mm_matIdx(test_name2)
            !print *,'mm_matIdx(test_name3) = ',mm_matIdx(test_name3)
            !print *,'mm_matIdx(test_name4) = ',mm_matIdx(test_name4)
            !print *,'mm_matIdx(test_name5) = ',mm_matIdx(test_name5)
            !print *,'mm_matIdx(test_name6) = ',mm_matIdx(test_name6)
            !print *,'enquiry_XS_test(1:6), 1.0/majorant_inv = ',enquiry_XS_test(1:6), 1.0/majorant_inv
            !print *,' '

              if (self % deform_type == 'swelling') then
               virtual_XS(1)=enquiry_XS_max / (self % y_factor * self % z_factor)
               virtual_XS(2)=enquiry_XS_max  / (self % x_factor * self % z_factor)
               virtual_XS(3)=enquiry_XS_max  / (self % x_factor * self % y_factor)
              elseif (self % deform_type == 'expansion') then
               virtual_XS(1)=enquiry_XS_max * self % z_factor
               virtual_XS(2)=enquiry_XS_max * self % z_factor
               virtual_XS(3)=enquiry_XS_max * self % z_factor
              end if
              virtual_XS_max_inv= ONE / maxval(virtual_XS)      ! Inverse of the maximum virtual cross-section!
              if (virtual_XS_max_inv < majorant_inv) then    ! Condition to check if the majorant will be exceeded in virtual perturbation!!
                 majorant_inv = virtual_XS_max_inv
              end if
          end if
    end if

    DTLoop:do
      distance = -log( p% pRNG % get() ) * majorant_inv
       if (self % Vir_Den_Swch==1) then
          if (trim(self % scale_type) == 'non_uniform') then
               if (any(self % pert_Mat_id == p % matIdx())) then
                self % x_factor_cur=self % x_factor
                self % y_factor_cur=self % y_factor
                self % z_factor_cur=self % z_factor
                self % deform_type = 'swelling'
               else
                 self % x_factor_cur=1.0
                 self % y_factor_cur=1.0
                 self % z_factor_cur=1.0
               end if
           end if

          if (((trim(self % scale_type) == 'non_uniform') .and. (any(self % pert_Mat_id == p % matIdx()))) &
              .or. (trim(self % scale_type) == 'uniform')) then    
                cosines(:) = p % dirGlobal()
                 !print *,'deform_type = ',self % deform_type
                 !print *,'direction_type = ',self % direction_type
                 !print *,'x_factor = ',self % x_factor_cur
                 !print *,'y_factor = ',self % y_factor_cur
                 !print *,'z_factor = ',self % z_factor_cur
                 !print *,'self % scale_type = ',self % scale_type
                 !print *,'pert_Mat_Id = ',self % pert_Mat_Id
                 !print *,'p % matIdx() = ',p % matIdx()
                 !print *,'OUTSIDE_FILL, VOID_MAT, UNDEF_MAT = ',OUTSIDE_FILL, VOID_MAT, UNDEF_MAT
                 !print *,''
                dist_vector=distance*cosines
               if (self % deform_type == 'swelling') then
                  Vir_Den_dist_vector(1)=dist_vector(1)*self % y_factor_cur*self % z_factor_cur
                  Vir_Den_dist_vector(2)=dist_vector(2)*self % x_factor_cur*self % z_factor_cur
                  Vir_Den_dist_vector(3)=dist_vector(3)*self % x_factor_cur*self % y_factor_cur
                  Vir_Den_dist=((Vir_Den_dist_vector(1)**2.0)+(Vir_Den_dist_vector(2)**(2.0))+(Vir_Den_dist_vector(3)**2.0))**0.5
                  flight_stretch_factor=Vir_Den_dist/distance
                  virtual_cosines(1)=cosines(1)*self % y_factor_cur*self % z_factor_cur/flight_stretch_factor
                  virtual_cosines(2)=cosines(2)*self % x_factor_cur*self % z_factor_cur/flight_stretch_factor
                  virtual_cosines(3)=cosines(3)*self % x_factor_cur*self % y_factor_cur/flight_stretch_factor
               elseif (self % deform_type == 'expansion') then  ! If the deformation type is Expansion!
                  Vir_Den_dist_vector(1)=dist_vector(1)/self % x_factor_cur
                  Vir_Den_dist_vector(2)=dist_vector(2)/self % y_factor_cur
                  Vir_Den_dist_vector(3)=dist_vector(3)/self % z_factor_cur
                  Vir_Den_dist=((Vir_Den_dist_vector(1)**2.0)+(Vir_Den_dist_vector(2)**(2.0))+(Vir_Den_dist_vector(3)**2.0))**0.5
                  flight_stretch_factor=Vir_Den_dist/distance
                  virtual_cosines(1)=cosines(1)/(self % x_factor_cur*flight_stretch_factor)
                  virtual_cosines(2)=cosines(2)/(self % y_factor_cur*flight_stretch_factor)
                  virtual_cosines(3)=cosines(3)/(self % z_factor_cur*flight_stretch_factor)
               else
                  print *,'Error in recognizing type of geometric deformation! Please check input!'
               end if
                call p % point(virtual_cosines)
                distance = Vir_Den_dist
                denizen_pert = 1    ! Indicates that this neutron has been found inside perturbed region! This condition remains unchanged until the loop ends!
               if (all(self % pert_Mat_id /= p % matIdx()) .and. (trim(self % scale_type) == 'non_uniform'))  then
                 print *,'deform_type = ',self % deform_type
                 print *,'direction_type = ',self % direction_type
                 print *,'scale_type = ',self % scale_type
                 print *,'dist_vector, dist = ',dist_vector,distance         ! Checked and values mathemematically feasible
                 print *,'Vir_Den_dist_vector = ',Vir_Den_dist_vector    ! Checked and values mathemematically feasible
                 print *,'Vir_Den_dist = ',Vir_Den_dist                  ! Checked and values mathemematically feasible 
                 print *,'flight_stretch_factor = ',flight_stretch_factor                ! Checked and values mathemematically feasible
                 print *,'virtual_cosines = ',virtual_cosines            ! Checked and values mathemematically feasible
                 print *,'x_factor_cur,y_factor_cur,z_factor_cur = ',self % x_factor_cur,self % y_factor_cur,self % z_factor_cur
                 print *,'pert_Mat_Id = ',self % pert_Mat_Id
                 print *,'p % matIdx() = ',p % matIdx()
                 print *,''
                end if
           end if
         if (trim(self % scale_type) == 'non_uniform') then
            if (any(self % pert_Mat_id == p % matIdx())) then
               pert_region_cond = 1         ! indicates neutron is in perturbed region
               mat_Idx_pert = p % matIdx()
               mat_Idx_unpert = 0
             else
               pert_region_cond = 0         ! indicates prisoner is in unperturbed region
               mat_Idx_unpert = p % matIdx()
               mat_Idx_pert = 0
            end if
         end if
       end if
      ! Move partice in the geometry

      call self % geom % teleport(p % coords, distance)

      if ((self % Vir_Den_Swch==1) .and. (trim(self % scale_type) == 'non_uniform')) then
         if ((pert_region_cond == 1) .and. (all(self % pert_Mat_id /= p % matIdx()))) then !&.and. (p % matIdx() /= OUTSIDE_FILL)) then
              pert_region_cond = 0
                !print *,'p % matIdx() = ',p % matIdx()
                !if (p % matIdx() == 0) then
                  !print *,'p % matIdx() = 0!'
                  !print *,'weight_pert_cond, p % w = ',weight_pert_cond, p % w
                !end if
                !print *,'weight_pert_cond outside = ',weight_pert_cond
              if ((abs(p % w - 1.0) < 1e-5) .or. (abs(p % w - 2.0) < 1e-5)) then ! prisoner is crossing over from its original perturbed cell to an unperturbed region!
                !print *,'denizen_pert, weight_pert_cond inside = ',denizen_pert, weight_pert_cond
                self % cross_over_cond = 1         ! markes the prisoner to have crossed over to a new region !
                p % w = p % w / (1.0 + ((self % z_factor - 1.0) / 1.0))   ! Weight is perturbed on cross-over to a new region
                !print *,'p % w = ',p % w
                !weight_pert_cond = 1
                !print *,'Mark 1 p % w = ',  p % w                                  ! Identifier that indicates this neutron's weight has been perturbed! This condition remains unchanged until the loop ends!
                !print *,'Prisoner just crossed over to unperturbed region, which is ',p % matIdx()
                !print *,'x_factor_cur,y_factor_cur,z_factor_cur = ',self % x_factor_cur,self % y_factor_cur,self % z_factor_cur
                !print *,'x_factor,y_factor,z_factor = ',self % x_factor,self % y_factor,self % z_factor
                !print *,'p % w,  mat_Idx_pert, p % matIdx() = ',p % w, mat_Idx_pert,p % matIdx()
                !print *,''
              elseif ((abs(p % w - (1.0 / (1.0 + ((self % z_factor - 1.0) / 1.0)))) < 1e-5) &
               .or. (abs(p % w - (2.0 / (1.0 + ((self % z_factor - 1.0) / 1.0)))) < 1e-5)) then ! prisoner has now gone back to his original region!
                self % cross_over_cond = 1         ! markes the prisoner to have crossed over to a new region !
                p % w = p % w * (1.0 + ((self % z_factor - 1.0) / 1.0))  ! Weight is perturbed on cross-over to a new region
                !weight_pert_cond = 0
                !print *,'Mark 2 p % w = ',p % w
              end if
                !print *,'Prisoner escaped to unperturbed region!'
                !print *,'mat_Idx_unpert,mat_Idx_pert, p % matIdx() = ',mat_Idx_unpert, mat_Idx_pert, p % matIdx()
                !print *,'p % w = ',p % w
                !print *,'pert_region_cond, self % cross_over_cond = ',pert_region_cond, self % cross_over_cond
                !print *,''
          elseif  ((pert_region_cond == 0) .and. (any(self % pert_Mat_id == p % matIdx()))) then
              pert_region_cond = 1
              if ((abs(p % w - 1.0) < 1e-5) .or. (abs(p % w - 2.0) < 1e-5)) then ! prisoner is crossing over from his original region to a new region !
                self % cross_over_cond = 1         ! markes the prisoner to have crossed over to a new region !
                p % w = p % w / (1.0 + ((self % z_factor - 1.0) / 1.0))   ! Weight is perturbed on cross-over to a new region
                !print *,'Prisoner just crossed over to perturbed region, which is ',p % matIdx()
                !print *,'x_factor_cur,y_factor_cur,z_factor_cur = ',self % x_factor_cur,self % y_factor_cur,self % z_factor_cur
                !print *,'p % w = ',p % w
                !print *,'x_factor,y_factor,z_factor = ',self % x_factor,self % y_factor,self % z_factor
                !print *,''
                !print *,'Mark 3 p % w = ',p % w
              elseif ((abs(p % w - (1.0 / (1.0 + ((self % z_factor - 1.0) / 1.0)))) < 1e-5) &
               .or. (abs(p % w - (2.0 / (1.0 + ((self % z_factor - 1.0) / 1.0)))) < 1e-5)) then ! prisoner has now gone back to his original region!
                self % cross_over_cond = 0         ! markes the prisoner is now back in his original region !
                p % w = p % w * (1.0 + ((self % z_factor - 1.0) / 1.0))    ! Weight perturbation is removed since the prisoner no longer contributes to the current
                !print *,' The prisoner is now back in perturbed region !'
                !print *,'Mark 4 p % w = ',p % w
             end if
                !print *,'Prisoner entered perturbed region!'
                !print *,'mat_Idx_unpert,mat_Idx_pert, p % matIdx() = ',mat_Idx_unpert, mat_Idx_pert, p % matIdx()
                !print *,'p % w = ',p % w
                !print *,'pert_region_cond, self % cross_over_cond = ',pert_region_cond, self % cross_over_cond
                !print *,''                    
          end if
            if ((p % w < 0.88) .and. ((self % z_factor - 1.0) > 1e-5)) then
                call fatalError(Here, "Warning! Incorrect current perturbation!")
                print *,'mat_Idx_unpert,mat_Idx_pert, p % matIdx() = ',mat_Idx_unpert, mat_Idx_pert, p % matIdx()
                print *,'p % w = ',p % w
                print *,'pert_region_cond, self % cross_over_cond = ',pert_region_cond, self % cross_over_cond
            end if                
       end if

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
            !if ((abs(p % w - 1.0 / (1.0 + ((self % z_factor - 1.0) / 2.0))) < 1e-5) .or. &
             !(abs(p % w - 2.0 / (1.0 + ((self % z_factor - 1.0) / 2.0))) < 1e-5)) then ! prisoner has now gone back to his original region!
                !p % w = p % w * (1.0 + ((self % z_factor - 1.0) / 2.0))
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

      !print *,'p % getCellIdx,sigmaT = ',p % getCellIdx()," ,",sigmaT
      !print *,'material_name = ',self % xsData % getTotalMatXS(p, mm_matIdx(p % matIdx()))
      ! Protect Against Sillines
      !if( sigmaT*majorant_inv < ZERO .or. ONE < sigmaT*majorant_inv) then
      !  call fatalError(Here, "TotalXS/MajorantXS is silly: "//numToChar(sigmaT*majorant_inv))
      !end if

    end do DTLoop
    
    call tally % reportTrans(p)
  end subroutine deltaTracking

  subroutine init(self, dict)
    class(transportOperatorDT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
         !Virtual Density Data call  begins !

      call init_super(self, dict)
      call dict % getorDefault(self % Vir_Den_Swch, 'vir_den_swch',0_shortInt)
       if (self % Vir_Den_Swch == 1) then
         call dict % getorDefault(self % deform_type, 'deform_type','swelling')
         call dict % getorDefault(self % direction_type, 'direction_type','isotropic')
         call dict % getorDefault(self % scale_type, 'scale','uniform')
         !print *,'Hello 1'
         !print *,'self % direction_type = ',self % direction_type
         !print *,'trim(self % direction_type) = ',trim(self % direction_type)

         if (trim(self % direction_type) == 'anisotropic') then
            print *,'Hello'
            call dict % getorDefault(self % x_factor, 'x_factor',1.0_defReal)
            call dict % getorDefault(self % y_factor, 'y_factor',1.0_defReal)
            call dict % getorDefault(self % z_factor, 'z_factor',1.0_defReal)
         else
            call dict % get(self % deform_factor, 'factor')
            self % x_factor=self % deform_factor
            self % y_factor=self % deform_factor
            self % z_factor=self % deform_factor
         end if

         self % x_factor_cur=self % x_factor
         self % y_factor_cur=self % y_factor
         self % z_factor_cur=self % z_factor
         print *,'Hello! Welcome to this virtual density perturbation simulation! Here are your input perturbation factors:'
         print *,'x_factor, y_factor, z_factor = ',self % x_factor, self % y_factor, self % z_factor
         print *,'Type of distortion = ', trim(self % direction_type), ', ', trim(self % deform_type)

         if (trim(self % scale_type) == 'non_uniform') then
               call dict % getorDefault(self % perturb_mat1, 'pert_mat_1','uniform')
               call dict % getorDefault(self % perturb_mat2, 'pert_mat_2','uniform')
               call dict % getorDefault(self % perturb_mat3, 'pert_mat_3','uniform')
               call dict % getorDefault(self % perturb_mat4, 'pert_mat_4','uniform')
               call dict % getorDefault(self % perturb_mat5, 'pert_mat_5','uniform')
               call dict % getorDefault(self % perturb_mat6, 'pert_mat_6','uniform')
               self % pert_Mat_Id(1) = mm_matIdx(self % perturb_mat1)
               self % pert_Mat_Id(2) = mm_matIdx(self % perturb_mat2)
               self % pert_Mat_Id(3) = mm_matIdx(self % perturb_mat3)
               self % pert_Mat_Id(4) = mm_matIdx(self % perturb_mat4)
               self % pert_Mat_Id(5) = mm_matIdx(self % perturb_mat5)
               self % pert_Mat_Id(6) = mm_matIdx(self % perturb_mat6)
               print *,'Perturbed MAT IDs are: ',self % pert_Mat_Id
         end if
        end if
!    Virtual Density Data call  ends !
  end subroutine init
  
end module transportOperatorDT_class
