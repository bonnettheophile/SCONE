module temporalPopClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState, P_PRECURSOR, P_NEUTRON
  use particleDungeon_class,      only : particleDungeon
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

  ! Nuclear Data Interface
  use nuclearDataReg_mod,     only : ndReg_get => get
  implicit none
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
  !! mytemporalPopClerk {
  !!   type temporalPopClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: temporalPopClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response

    ! Useful data
    integer(shortInt)  :: width = 0

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportTemporalPopIn
    procedure :: reportTemporalPopOut
    procedure :: reportCycleStart

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type temporalPopClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(temporalPopClerk), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i

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

    ! Get names of response dictionaries
    call dict % get(responseNames,'response')

    ! Load responses
    allocate(self % response(size(responseNames)))
    do i=1, size(responseNames)
      call self % response(i) % init(dict % getDictPtr( responseNames(i) ))
    end do

    ! Set width
    self % width = size(responseNames)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(temporalPopClerk), intent(inout) :: self

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
    class(temporalPopClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleStart_CODE, temporalPop_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(temporalPopClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = size(self % response)
    if(allocated(self % map)) S = S * self % map % bins(0)

  end function getSize

  subroutine reportCycleStart(self, start, mem)
    class(temporalPopClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: start
    type(scoreMemory), intent(inout)      :: mem
    class(nuclearDatabase), pointer        :: xsData
    integer(shortInt)                      :: binIdx, i, j
    integer(longInt)                       :: adrr
    real(defReal)                          :: scoreVal
    type(particleState)                    :: state
    type(particle)                         :: p
    character(100), parameter :: Here =' reportCycleStart (temporalPopClerk_class.f90)'

    ! Append all bins
    do i=1, self % width

      do j = 1, start % popSize()
        ! Get current particle state
        call start % copy(p, j)

        if (.not. p % isdead .and. p % type == P_NEUTRON) then
          state = p
          state % timeBinIdx = p % timeBinIdx + 1
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
          scoreVal = p % w * self % response(i) % get(p, xsData)
          call mem % score(scoreVal, adrr + i)

        end if
        
      end do
      
    end do


  end subroutine reportCycleStart

  !!
  !! Process temporal population report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTemporalPopIn(self, p, mem)
    class(temporalPopClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    type(scoreMemory), intent(inout)       :: mem
    class(nuclearDatabase), pointer        :: xsData
    integer(shortInt)                      :: binIdx, i
    integer(longInt)                       :: adrr
    real(defReal)                          :: scoreVal
    type(particleState)                    :: state
    character(100), parameter :: Here =' reportTemporalPopIn (temporalPopClerk_class.f90)'

    if (p % type == P_PRECURSOR) then
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
        scoreVal = p % w * self % response(i) % get(p, xsData)
        call mem % score(scoreVal, adrr + i)

      end do
    end if

  end subroutine reportTemporalPopIn

  !!
  !! Process temporal population report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTemporalPopOut(self, p, mem)
    class(temporalPopClerk), intent(inout) :: self
    class(particle), intent(in)            :: p
    type(scoreMemory), intent(inout)       :: mem
    class(nuclearDatabase), pointer        :: xsData
    integer(shortInt)                      :: binIdx, i
    integer(longInt)                       :: adrr
    real(defReal)                          :: scoreVal
    type(particleState)                    :: state
    character(100), parameter :: Here =' reportTemporalPopOut (temporalPopClerk_class.f90)'

    if (p % type == P_PRECURSOR) then
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
        scoreVal = -p % w * self % response(i) % get(p, xsData)
        call mem % score(scoreVal, adrr + i)

      end do
    end if

  end subroutine reportTemporalPopOut

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(temporalPopClerk), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'temporalPopClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(temporalPopClerk), intent(in)        :: self
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
    if(allocated(self % map)) then
      resArrayShape = [size(self % response), self % map % binArrayShape()]
    else
      resArrayShape = [size(self % response)]
    end if

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i=1,product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val, std)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

end module temporalPopClerk_class
