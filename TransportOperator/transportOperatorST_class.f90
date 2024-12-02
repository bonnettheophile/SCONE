!!
!! Transport operator for surface tracking
!!
module transportOperatorST_class
  use numPrecision
  use universalVariables

  use errors_mod,               only : fatalError
  use particle_class,           only : particle
  use particleDungeon_class,    only : particleDungeon
  use dictionary_class,         only : dictionary

  ! Superclass
  use transportOperator_inter,  only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,           only : geometry, distCache

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,         only : tallyAdmin

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
    logical(defBool)                 :: virtual_density, cross_over = .false.
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
    integer(shortInt)                         :: event, collFate, i, current_mat
    real(defReal)                             :: sigmaT, dist, sigmaTrack, invSigmaTrack, &
                                                 speed, time, virtual_dist, flight_stretch_factor
    real(defReal), dimension(3)                :: preCosines, cosines, virtual_cosines, real_vector, virtual_vector
    type(distCache)                           :: cache
    real(defReal), parameter                  :: tol  = 1.0E-12
    character(100), parameter:: Here = 'surfaceTracking (transportOperatorST_class.f90)'

    preCosines = p % dirGlobal()

    
    STLoop: do
        
      ! Get local conditions
      call self % localConditions(p)
      
      sigmaTrack = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)

      ! Obtain the local cross-section, depending on the material
      ! This branch is called in the case of voids with no imposed XS
      if (sigmaTrack < tol) then
        
        dist = INFINITY
        invSigmaTrack = INFINITY
        sigmaT = ZERO

      else
        
        invSigmaTrack = ONE / sigmaTrack
        dist = -log( p % pRNG % get()) * invSigmaTrack
      
        ! Obtain the local cross-section
        sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())

        ! Should never happen  ! Catches NaN distances
        if (dist /= dist) then
          print *, "Particle location: ", p % rGlobal()
          print *, "Particle direction: ", p % dirGlobal()
          print *, "Total XS: ", sigmaT
          call fatalError(Here, "Distance is NaN")
        end if

      end if

      speed = p % getSpeed()
      time = dist / speed + p % time
      
      ! Set a max flight distance due to hitting the time-boundary
      if (p % timeMax > ZERO .and. time > p % timeMax) then
        dist = speed * (p % timeMax - p % time)
        collFate = AGED_FATE
      else
        collFate = NO_FATE
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
          cosines(:) = preCosines 
          real_vector = dist*cosines

          if (self % deform_type(current_mat) == 'swelling') then
            virtual_vector(1) = real_vector(1) * p % f(2) * p % f(3)
            virtual_vector(2) = real_vector(2) * p % f(1) * p % f(3)
            virtual_vector(3) = real_vector(3) * p % f(1) * p % f(2)
            virtual_dist = sqrt(sum(virtual_vector**2))
            flight_stretch_factor = virtual_dist/dist
            virtual_cosines(1) = cosines(1) * p % f(2) * &
                p % f(3) / flight_stretch_factor
            virtual_cosines(2) = cosines(2) * p % f(1) * &
                p % f(3) / flight_stretch_factor
            virtual_cosines(3) = cosines(3) * p % f(1) * &
                p % f(2) / flight_stretch_factor
          elseif (self % deform_type(current_mat) == 'expansion') then
            virtual_vector = real_vector/p % f
            virtual_dist = sqrt(sum(virtual_vector**2))
            flight_stretch_factor = virtual_dist/dist
            virtual_cosines = cosines / (p % f*flight_stretch_factor)
          else
            call fatalError(Here, 'Unrecognised geometric deformation')
          end if
        
          call p % point(virtual_cosines)
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

      ! Advance in time
      p % time = p % time + dist / speed

      ! Set fate if a collision occurred
      if (event == COLL_EV) p % fate = collFate
      
      call p % point(preCosines)
      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      select case(p % matIdx())
      
        ! Kill particle if it has leaked
        case(OUTSIDE_FILL)
          p % isDead = .true.
          p % fate = LEAK_FATE

        ! Give error if the particle somehow ended in an undefined material
        case(UNDEF_MAT)
          print *, "Particle location: ", p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")
        
        ! Give error if the particle is in a region with overlapping cells
        case(OVERLAP_MAT)
          print *, "Particle location: ", p % rGlobal()
          call fatalError(Here, "Particle is in overlapping cells")
        
        case default
          ! All is well

      end select

      ! Return if particle is stopped by death, or aging
      if (p % isDead .or. p % fate == AGED_FATE) exit STLoop

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (event == COLL_EV) then
        if (p % pRNG % get() < sigmaT*invSigmaTrack) then
          exit STLoop
        else
          call tally % reportInColl(p, .true.)
        end if
      end if

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

    ! Initialise cache
    if (dict % isPresent('cache')) then
      call dict % get(self % cache, 'cache')
    end if

  end subroutine init

end module transportOperatorST_class
