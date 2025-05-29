module importanceClerk_class
    use numPrecision
    use tallyCodes
    use endfConstants
    use universalVariables
    use genericProcedures,          only : fatalError
    use dictionary_class,           only : dictionary
    use particle_class,             only : particle, particleState, ancestorState
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

    use particleDungeon_class,      only : particleDungeon

    implicit none
    private

    !! 
    !! Importance estimator
    !! Calculate neutron importance through Serpent 2 formula (see paper 2014)
    !! Assumes the dungeon to be sorted in ascending order of oldest ancestor
    !!
    !! Private Members:
        !!   filter   -> Space to store tally Filter
        !!   map      -> Space to store tally Map
    !!
    !! Interface
    !!   tallyClerk Interface
    !! SAMPLE DICTIONATY INPUT:
    !!
    !! myImportanceClerk {
    !!   type importanceClerk;
    !!   # filter { <tallyFilter definition> } #
    !!   # map    { <tallyMap definition>    } #
    !! }
    !!

type, public, extends(tallyClerk) :: importanceClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type importanceClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(importanceClerk), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name

    ! Assign name
    call self % setName(name)

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(importanceClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

  end subroutine kill

    !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(importanceClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(importanceClerk), intent(in) :: self
    integer(shortInt)                 :: S

    if(allocated(self % map)) S = self % map % bins(0)

  end function getSize

  !!
  !! Process beginning of a cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(importanceClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)      :: end
    type(scoreMemory), intent(inout)        :: mem
    type(particleState)                     :: state, ancestorState
    integer(shortInt)                       :: idx, lastID, binIdx
    real(defReal)                           :: scoreVal
    integer(longInt)                        :: addr

    ! Get first particle for testing if enough generation passed for scoring importance
    state = end % get(1)
    scoreVal = 0
    if (all(state % ifpAncestorsID /= 0)) then
      ! Initialise lastID for detecting change of oldest ancestor
      lastID = state % ifpAncestorsID(1)
      do idx = 1, end % popSize()
        ! Get current particle state
        state = end % get(idx)
        ! Get ancestor state to determine mapping idx
        ancestorState = end % getAncestor(idx, 1)

        ! Check if ancestor is within filter
        if (allocated(self % filter)) then
            if (self % filter % isFail(ancestorState)) return
        end if

        ! Find bin index
        if (allocated(self % map)) then
            binIdx = self % map % map(ancestorState)
        else
            binIdx = 1
        end if

        ! Return if invalid bin index
        if (binIdx == 0) return

        ! Add score if same oldest ancestor as previous particle
        if (state % ifpAncestorsID(1) == lastID) then
          scoreVal = scoreVal + ONE * state % ifpAncestorWgt(1)
        else 
          ! Score if different oldest ancestor from previous particle
          addr = self % getMemAddress() + (binIdx - 1)  - 1
          call mem % score(scoreVal, addr)
          ! Start scoring for next oldest ancestor
          scoreVal = ONE * state % ifpAncestorWgt(1)
          ! Update lastId
          lastID = state % ifpAncestorsID(1)
        end if
      end do
    end if
  end subroutine reportCycleEnd

!!
!! Display convergence progress on the console
!!
!! See tallyClerk_inter for details
!!
  subroutine display(self, mem)
  class(importanceClerk), intent(in)  :: self
  type(scoreMemory), intent(in)      :: mem

    print *, 'importanceClerk does not support display yet'

  end subroutine display

    !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(importanceClerk), intent(in)          :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If importance clerk has map print map information
    if (allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    if (allocated(self % map)) then
      resArrayShape = [self % map % binArrayShape()]
    else
      resArrayShape = [1]
    end if

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val,std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

end module importanceClerk_class