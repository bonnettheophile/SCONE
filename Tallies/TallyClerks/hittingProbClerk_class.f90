module hittingProbClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState, P_NEUTRON
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot

  ! Particle Dungeon
  use particleDungeon_class,  only : particleDungeon

  ! Nuclear Data Interface
  use nuclearDataReg_mod,     only : ndReg_get => get
  implicit none

  integer(shortInt), parameter, public :: UPPER = 1, LOWER = 2

  private

  !!
  !! temporal estimator of populations
  !! Calculates discrete populations at time boundaries
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!   response -> Array of responses
  !!   width    -> Number of responses (# of result bins for each map position)
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myhittingProbClerk {
  !!   type hittingProbClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: hittingProbClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response
    real(defReal)                                    :: maxT
    integer(shortInt)                                :: maxPop, minPop
    logical(defBool), dimension(:), allocatable      :: firstHitsNeutron
    real(defReal), dimension(:), allocatable         :: hittingTimesNeutron
    real(defReal), dimension(:), allocatable         :: currentNeutronPops
    integer(shortInt)                                :: cycles
    integer(shortInt)                                :: thresholdKind
    type(particleDungeon), allocatable               :: neutrons

    ! Useful data
    integer(shortInt)  :: width = 0

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportHittingProbIn
    procedure  :: reportHittingProbOut
    procedure  :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print
  end type hittingProbClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(hittingProbClerk), intent(inout)      :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i, j
    class(dictionary), pointer                  :: responseDict
    character(nameLen)                          :: type
    character(100), parameter :: Here = 'init (hittingProbClerk_class.f90)'

    ! Assign name
    call self % setName(name)

    ! Load filter
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

    ! Get max time to consider
    call dict % get(self % maxT,'maxT')

    ! Get if its a lower bound or upper bound
    call dict % getOrDefault(type, 'kind', 'upper')
    ! Get population max threshold
    if (trim(type) == 'upper') then
      self % thresholdKind = UPPER
      call dict % get(self % maxPop,'maxPop')
    ! Get population min threshold
    else if (trim(type) == 'lower') then
      self % thresholdKind = LOWER
      call dict % get(self % minPop, 'minPop')
    else
      call fatalError(Here, 'Invalid input for threshold kind')
    end if
    



    ! Get population threshold
    call dict % get(self % cycles,'cycles')

    ! Get names of response dictionaries
    call dict % get(responseNames,'response')

    ! Load responses
    allocate(self % response(size(responseNames)))
    do i=1, size(responseNames)
      call self % response(i) % init(dict % getDictPtr( responseNames(i) ))

      responseDict => dict % getDictPtr( responseNames(i) )
      call responseDict % get(type,'type')

      if (type == "neutronResponse") then
        allocate(self % neutrons)
        call self % neutrons % init(20)
        allocate(self % currentNeutronPops(self % cycles))
        allocate(self % firstHitsNeutron(self % cycles))
        allocate(self % hittingTimesNeutron(self % cycles))
        do j=1, self % cycles
          self % currentNeutronPops(j) = ZERO
          self % firstHitsNeutron(j) = .true.
          self % hittingTimesNeutron(j) = -ONE
        end do
      end if
    end do

    ! Set width
    self % width = size(responseNames)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(hittingProbClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if(allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if(allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    ! Kill and deallocate responses
    if(allocated(self % response)) then
      deallocate(self % response)
    end if

    self % width = 0

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(hittingProbClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [hittingProb_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(hittingProbClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = size(self % response)
    if(allocated(self % map)) S = S * self % map % bins(0)

  end function getSize

  !!
  !! Process temporal population report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHittingProbIn(self, p, mem)
    class(hittingProbClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    type(scoreMemory), intent(inout)       :: mem
    class(nuclearDatabase), pointer        :: xsData
    integer(shortInt)                      :: i, batchIdx, binIdx
    real(defReal)                          :: scoreVal
    type(particleState)                    :: state
    integer(longInt)                       :: adrr
    character(100), parameter :: Here =' reporthittingProbIn (hittingProbClerk_class.f90)'

    batchIdx = mod(mem % batchN, self % cycles) + 1
    if (self % firstHitsNeutron(batchIdx) .eqv. .false.) return

    if (p % time <= self % maxT) then
      ! Append all bins
      do i=1, self % width

        ! Get current particle state
        state = p

        ! Check if within filter
        if(allocated( self % filter)) then
          if(self % filter % isFail(state)) return
        end if

        ! Find bin index
        if(allocated(self % map)) then
          binIdx = self % map % map(state)
        else
          binIdx = 1
        end if

        ! Return if invalid bin index
        if (binIdx == 0) return

        ! Calculate bin address
        adrr = self % getMemAddress() + self % width * (binIdx - 1) - 1

        xsData => ndReg_get(p % getType(), where = Here)
        scoreVal = p % w * self % response(i) % get(p, xsData) !self % response(i) % get(p, xsData)

        if (scoreVal == ZERO) return
        if (allocated(self % neutrons) .and. p % type == P_NEUTRON) then
          call mem % score(scoreVal, adrr + i)
        end if

      end do
    end if

  end subroutine reportHittingProbIn

  !!
  !! Process temporal population report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHittingProbOut(self, p, mem)
    class(hittingProbClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    type(scoreMemory), intent(inout)       :: mem
    class(nuclearDatabase), pointer        :: xsData
    integer(shortInt)                      :: i, batchIdx, binIdx
    real(defReal)                          :: scoreVal
    type(particleState)                    :: state
    integer(longInt)                       :: adrr
    character(100), parameter :: Here =' reporthittingProbOut (hittingProbClerk_class.f90)'

    batchIdx = mod(mem % batchN, self % cycles) + 1
    if (self % firstHitsNeutron(batchIdx) .eqv. .false.) return

    if (p % time <= self % maxT) then
      ! Append all bins
      do i=1, self % width

        ! Get current particle state
        state = p

        ! Check if within filter
        if(allocated( self % filter)) then
          if(self % filter % isFail(state)) return
        end if

        ! Find bin index
        if(allocated(self % map)) then
          binIdx = self % map % map(state)
        else
          binIdx = 1
        end if

        ! Return if invalid bin index
        if (binIdx == 0) return

        ! Calculate bin address
        adrr = self % getMemAddress() + self % width * (binIdx - 1) - 1

        xsData => ndReg_get(p % getType(), where = Here)
        scoreVal = - p % w * self % response(i) % get(p, xsData)
        if (scoreVal == ZERO) return
        if (allocated(self % neutrons) .and. p % type == P_NEUTRON) then
          call mem % score(scoreVal, adrr + i)
        end if

      end do
    end if

  end subroutine reportHittingProbOut

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(hittingProbClerk), intent(inout)  :: self
    class(particleDungeon), intent(in)      :: end
    type(scoreMemory), intent(inout)        :: mem
    integer(shortInt)                       :: batchIdx, i
    integer(longInt)                        :: adrr
    real(defReal)                           :: accScore
    integer(shortInt),dimension(:),allocatable :: resArrayShape

    batchIdx = mod(mem % batchN, self % cycles) + 1

    if(allocated(self % map)) then
      resArrayShape = [size(self % response), self % map % binArrayShape()]
    else
      resArrayShape = [size(self % response)]
    end if

    do i = 1, product(resArrayShape)
      adrr = self % getMemAddress() + self % width * (i - 1) - 1
      accScore = mem % getScore(adrr)
      self % currentNeutronPops(batchIdx) = accScore
      if (self % currentNeutronPops(batchIdx) >= self % maxPop .and. self % thresholdKind == UPPER) then
        self % firstHitsNeutron(batchIdx) = .false.
        self % hittingTimesNeutron(batchIdx) = i
      else if (self % currentNeutronPops(batchIdx) <= self % minPop .and. self % thresholdKind == LOWER) then
        self % firstHitsNeutron(batchIdx) = .false.
        self % hittingTimesNeutron(batchIdx) = i
      end if

    end do

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(hittingProbClerk), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'hittingProbClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(hittingProbClerk), intent(in)        :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If temporal population clerk has map print map information
    if( allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array

    resArrayShape = [size(self % response)]

    ! Start array
    name ='HittingProbability'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    val = ZERO
    do i=1, self % cycles
      if (self % hittingTimesNeutron(i) >= ZERO) val = val + ONE
    end do
    val = val / self % cycles
    call outFile % addValue(val)
    call outFile % endArray()

    ! Start array
    resArrayShape = [self % cycles]
    name ='HittingTimeIdxs'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i=1, product(resArrayShape)
      val = self % hittingTimesNeutron(i)
      call outFile % addValue(val)
    end do

    call outFile % endArray()

    ! Write results.
    ! Get shape of result array
    if(allocated(self % map)) then
      resArrayShape = [size(self % response), self % map % binArrayShape()]
    else
      resArrayShape = [size(self % response)]
    end if

    ! Start array
    name ='pop'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i=1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val, std)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

end module hittingProbClerk_class
