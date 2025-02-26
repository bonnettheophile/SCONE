!!
!! Transport operator for surface tracking
!!
module transportOperatorST_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use coord_class,                only : coordList, coord

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry, distCache

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
  !! Transport operator that moves a particle with surface tracking
  !!
  !! Sample Input Dictionary:
  !!   trans { type transportOperatorST; cache 0; }
  !!
  type, public, extends(transportOperator):: transportOperatorST
    logical(defBool)  :: cache = .true.
    real(defReal)                    :: product_factor 
    real(defReal), dimension(:,:), allocatable      :: vector_factor
    logical(defBool)                 :: virtual_density, cross_over = .false., swelled
    character(nameLen)               :: direction_type, scale_type
    character(nameLen), allocatable  :: pert_mat(:), deform_type(:)
    integer(shortInt), allocatable   :: pert_mat_id(:)
    integer(shortInt)                :: nb_pert_mat
  contains
    procedure:: transit => surfaceTracking
    ! Override procedure
    procedure:: init
  end type transportOperatorST

contains

  !!
  !! Performs surface tracking until a collision point is found
  !!
  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorST), intent(inout):: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)      :: thisCycle
    class(particleDungeon), intent(inout)      :: nextCycle
    integer(shortInt)                         :: event, i, current_mat
    real(defReal)                             :: sigmaT, dist, virtual_dist, flight_stretch_factor
    real(defReal), dimension(3)                :: cosines, virtual_cosines, real_vector, virtual_vector
    type(distCache)                           :: cache
    logical(defBool)                          :: repoint, first_flight
    character(100), parameter:: Here = 'surfaceTracking (transportOperatorST_class.f90)'
    first_flight = .true.
    STLoop: do

      ! Obtain the local cross-section
      if (p % matIdx() == VOID_MAT) then
        dist = INFINITY

      else
        sigmaT = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)
        dist = -log( p % pRNG % get()) / sigmaT*self % product_factor

        ! Should never happen  ! Catches NaN distances
        if (dist /= dist) then 
          print *, p % matIdx()
          call fatalError(Here, "Distance is NaN")
        end if

      end if

      if (self % virtual_density) then
        ! If uniform virtual density, should always be 1
        current_mat = 1
        if (trim(self % scale_type) == 'non_uniform') then
          p % isPerturbed = .false.
          do i = 1, self % nb_pert_mat 
            if (self % pert_mat_id(i) == p % matIdx()) then
              p % isPerturbed = .true. ! Set particle to be perturbed
              current_mat = i  ! Set current perturbated material
              exit
            end if
          end do
        end if

        if (p % isPerturbed .or. trim(self % scale_type) == 'uniform') then
          cosines(:) = p % dirGlobal()
          real_vector = dist*cosines
          if (self % deform_type(current_mat) == 'swelling') then
            virtual_vector(1) = real_vector(1) * self % vector_factor(2, current_mat) * self % vector_factor(3, current_mat)
            virtual_vector(2) = real_vector(2) * self % vector_factor(1, current_mat) * self % vector_factor(3, current_mat)
            virtual_vector(3) = real_vector(3) * self % vector_factor(1, current_mat) * self % vector_factor(2, current_mat)
            virtual_dist = sqrt(sum(virtual_vector**2))
            flight_stretch_factor = virtual_dist/dist
            virtual_cosines(1) = cosines(1) * self % vector_factor(2, current_mat) * &
                self % vector_factor(3, current_mat) / flight_stretch_factor
            virtual_cosines(2) = cosines(2) * self % vector_factor(1, current_mat) * &
                self % vector_factor(3, current_mat) / flight_stretch_factor
            virtual_cosines(3) = cosines(3) * self % vector_factor(1, current_mat) * &
                self % vector_factor(2, current_mat) / flight_stretch_factor
          elseif (self % deform_type(current_mat) == 'expansion') then
            virtual_vector = real_vector/self % vector_factor(:,current_mat) 
            virtual_dist = sqrt(sum(virtual_vector**2))
            flight_stretch_factor = virtual_dist/dist
            virtual_cosines = cosines / (self % vector_factor(:,current_mat)*flight_stretch_factor)
          else
            call fatalError(Here, 'Unrecognised geometric deformation')
          end if
        
          if (first_flight .and. trim(self % scale_type) == 'uniform') then
            call p % point(virtual_cosines)
            first_flight = .false.
          else if (first_flight .and. trim(self % scale_type) == 'non_uniform') then
            call p % point(virtual_cosines)
            first_flight = .false.
          end if
          dist = virtual_dist
          
        end if
      end if

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop.
      if (self % cache) then
        call self % geom % move_withCache(p % coords, dist, event, cache)

      else
        call self % geom % move(p % coords, dist, event)

      end if
      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      if (self % virtual_density .and. trim(self % scale_type) == 'non_uniform') then 
        p % lastPerturbed = p % isPerturbed
        repoint = .false.
        if (any(self % pert_mat_id == p % matIdx())) then
          p % isPerturbed = .true.
        else
          p % isPerturbed = .false.
        end if

        ! If crossing to unperturbed region, recover non perturbed direction
        if ( (p % lastPerturbed .and. (.not. p % isPerturbed))) then 
          call p % point(cosines)
          first_flight = .true.
        end if
      end if


      ! Kill particle if it has leaked
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = LEAK_FATE
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      ! Return if particle stoped at collision (not cell boundary)
      if (event == COLL_EV .or. p % isDead) exit STLoop

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking

  !!
  !! Initialise ST operator from a dictionary
  !!
  !! See transportOperator_inter for details
  !!
  subroutine init(self, dict)
    class(transportOperatorST), intent(inout):: self
    class(dictionary), intent(in)             :: dict
    character(nameLen)                        :: input
    real(defReal), allocatable, dimension(:)  :: vec                           
    integer(shortInt)                         :: index

    ! Initialise superclass
    call init_super(self, dict)

    ! Initialise virtual density
    call dict % getorDefault(self % virtual_density, 'virtual_density', .false.)
    call dict % getorDefault(self % swelled, 'swelled', .false.)
    self % product_factor = ONE
    if (self % swelled) then
        allocate(self % deform_type(1))
        allocate(self % vector_factor(3, 1))
        call dict % get(self % deform_type(1), "deform_type_1")
        call dict % get(vec, "factor_1")
        self % vector_factor(:,1) = vec
        self % product_factor = product(self % vector_factor)
    end if
    if (self % virtual_density) then
      call dict % getorDefault(self % direction_type, 'direction_type','isotropic')
      call dict % getorDefault(self % scale_type, 'scale','uniform')
      call dict % getorDefault(self % nb_pert_mat, 'nb_pert_mat', 1)

      !self % product_factor = self % vector_factor(1) * self % vector_factor(:,2) * self % vector_factor(3)

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
          print *, self % pert_mat_id(index), self % pert_mat(index)

          input = 'deform_type_'
          write(input, '(I0)') index
          input = trim('deform_type_')//trim(input)
          call dict % get(self % deform_type(index), trim(input))
          print *, self % deform_type(index)
        end do
      else
        allocate(self % deform_type(1))
        allocate(self % vector_factor(3, 1))
        call dict % get(self % deform_type(1), "deform_type_1")
        call dict % get(vec, "factor_1")
        self % vector_factor(:,1) = vec
      end if
    end if

    ! Initialise cache
    if (dict % isPresent('cache')) then
      call dict % get(self % cache, 'cache')
    end if

  end subroutine init

end module transportOperatorST_class
