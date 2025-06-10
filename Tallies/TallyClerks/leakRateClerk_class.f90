module leakRateClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use endfConstants
    use genericProcedures,          only : fatalError
    use dictionary_class,           only : dictionary
    use particle_class,             only : particle, particleState
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
  
    implicit none
    private
  
    !!
    !! Analog estimator of rate of leakage
    !! Calculates number of leaking particles
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
    !! myLeakageClerk {
    !!   type leakRateClerk;
    !!   # filter { <tallyFilter definition> } #
    !!   # map    { <tallyMap definition>    } #
    !! }
    !!
    type, public, extends(tallyClerk) :: leakRateClerk
      private
      ! Filter, Map & Vector of Responses
      class(tallyFilter), allocatable                  :: filter
      class(tallyMap), allocatable                     :: map
  
      ! Useful data
      logical(defBool)   :: normByPop = .false.
      real(defReal)      :: invPopSize = ONE ! 1 so that without normalization it behaves as usual
      integer(shortInt)  :: width = 0
  
    contains
      ! Procedures used during build
      procedure  :: init
      procedure  :: kill
      procedure  :: validReports
      procedure  :: getSize
  
      ! File reports and check status -> run-time procedures
      procedure  :: reportHist
  
      ! Output procedures
      procedure  :: display
      procedure  :: print
  
    end type leakRateClerk
  
  contains
  
    !!
    !! Initialise clerk from dictionary and name
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine init(self, dict, name)
      class(leakRateClerk), intent(inout)        :: self
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
  
      self % width = 1

    end subroutine init
  
    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(leakRateClerk), intent(inout) :: self
  
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
  
      self % width   = 0

    end subroutine kill
  
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    !! See tallyClerk_inter for details
    !!
    function validReports(self) result(validCodes)
      class(leakRateClerk),intent(in)           :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
  
      validCodes = [hist_CODE]
  
    end function validReports
  
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!
    elemental function getSize(self) result(S)
      class(leakRateClerk), intent(in) :: self
      integer(shortInt)                :: S
  
      if (allocated(self % map)) S = self % width * self % map % bins(0)
  
    end function getSize
  
    !!
    !!
    !!
  
    subroutine reportHist(self, p, xsData, mem)
        class(leakRateClerk), intent(inout) :: self
        class(particle), intent(in)             :: p
        class(nuclearDatabase),intent(inout)    :: xsData
        type(scoreMemory), intent(inout)        :: mem
        real(defReal)                           :: histWgt
        integer(longInt)                        :: adrr
        type(particleState)                    :: state
        integer(shortInt)                       :: binIdx
    
        if( p % fate == leak_FATE) then
          ! Obtain and score history weight
          histWgt = p % w
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
          adrr = self % getMemAddress() + self % width * (binIdx -1)  - 1
          ! Score analog leakage
          call mem % score( histWgt, adrr)
    
        end if
    
      end subroutine reportHist

    !!
    !! Display convergance progress on the console
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine display(self, mem)
      class(leakRateClerk), intent(in)  :: self
      type(scoreMemory), intent(in)      :: mem
  
      print *, 'leakRateClerk does not support display yet'
  
    end subroutine display
  
    !!
    !! Write contents of the clerk to output file
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine print(self, outFile, mem)
      class(leakRateClerk), intent(in)          :: self
      class(outputFile), intent(inout)           :: outFile
      type(scoreMemory), intent(in)              :: mem
      real(defReal)                              :: val, std
      integer(shortInt)                          :: i
      integer(shortInt),dimension(:),allocatable :: resArrayShape
      character(nameLen)                         :: name
  
      ! Begin block
      call outFile % startBlock(self % getName())
  
      ! If collision clerk has map print map information
      if( allocated(self % map)) then
        call self % map % print(outFile)
      end if
  
      ! Write results.
      ! Get shape of result array
      if(allocated(self % map)) then
        resArrayShape = [1, self % map % binArrayShape()]
      else
        resArrayShape = [1]
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
  
  end module leakRateClerk_class
  