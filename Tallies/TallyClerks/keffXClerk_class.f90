module keffXClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use endfConstants
    use dictionary_class,      only : dictionary
    use genericProcedures,     only : fatalError
    use particle_class,        only : particle, particleState
    use particleDungeon_class, only : particleDungeon
    use outputFile_class,      only : outputFile
    use legendrePoly_func,     only : evaluateLegendre
  
  ! Nuclear Data Interfaces
    use nuclearDataReg_mod,         only : ndReg_get => get
    use nuclearDatabase_inter,      only : nuclearDatabase
    use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast
    use neutronXSPackages_class,    only : neutronMacroXSs

    use scoreMemory_class,     only : scoreMemory
    use tallyResult_class,     only : tallyResult, tallyResultEmpty
    use tallyClerk_inter,      only : tallyClerk, kill_super => kill

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap
  
    implicit none
    private
    
    !! 
    !! Estimator for polynomial chaos coefficients when using virtual densities
    !! Output is a 2D array of size (order+1) ** number. We start assuming a single random parameter
    !! Additional assumption: population control is enforcing starting weight is constant

    !! Private members:
    !!   chaosOfPop   -> Array for end of generation chaos model
    !!   map          -> Space to store tally map
    !!   P            -> Order of chaos model
    !!   startPop     -> Starting population for each cycle

    !! Interface
    !!   tallyClerk Interface

    !! SAMPLE DICTIONARY INPUT
    !! 
    !! myCoeffofChaosClerk {
    !!    type keffXClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: keffXClerk
      private
      class(tallyMap), allocatable                     :: map
      integer(shortInt)                                :: width = 0

      
    contains
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize

      procedure :: reportHist
      procedure :: reportCycleEnd

      procedure :: print
      procedure :: display

    end type keffXClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(keffXClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      character(100),parameter :: Here = 'init (keffXClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if
    self % width = 3
    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(keffXClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      ! Kill and deallocate map
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      
      

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(keffXClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [hist_CODE, cycleEnd_CODE]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(keffXClerk), intent(in) :: self
      integer(shortInt)                    :: S

      S = self % width
      if(allocated(self % map)) S = 3 * self % map % bins(0)
    end function getSize


    subroutine reportCycleEnd(self, end, mem)
      class(keffXClerk), intent(inout) :: self
      class(particleDungeon), intent(in)    :: end
      type(scoreMemory), intent(inout)      :: mem
      integer(shortInt)                     :: N, i
      integer(longInt)                      :: addr
      real(defReal)                         :: value, events, score

      N = self % map % bins(0)
      do i = 1, N
        addr = self % getMemAddress() + self % width * (i - 1) - 1
        value = mem % getScore(addr + 1)
        events = mem % getScore(addr + 2)

        if (events /= ZERO) then 
          score = value / events
        else
          score = ZERO
        end if

        call mem % score(score, addr + 3)
        !call mem % setScore(addr + 1, ZERO)
        !call mem % setScore(addr + 2, ZERO)
      end do

      
    end subroutine reportCycleEnd
    !!
    !! Process history report
    !! Gets fate code from the particle
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportHist(self, p, xsData, mem)
      class(keffXClerk), intent(inout) :: self
      class(particle), intent(in)             :: p
      class(nuclearDatabase),intent(inout)    :: xsData
      type(scoreMemory), intent(inout)        :: mem
      real(defReal)                           :: score
      integer(shortInt)                       :: binIdx
      integer(longInt)                        :: addr
      type(particleState)                     :: state


      if (p % fate /= leak_FATE) then
        score = p % w
      else
        score = ZERO
      end if
      state = p

      if (allocated(self % map)) then
        binIdx = self % map % map(state)
      else
        binIdx = 1
      end if
      
      ! Return if invalid bin index
      if (binIdx == 0) return

      addr = self % getMemAddress() + self % width * (binIdx - 1)  - 1

      call mem % score(score, addr + 1)
      call mem % score(ONE, addr + 2)
    end subroutine reportHist

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(keffXClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'keffXClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(keffXClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: i
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        character(nameLen)                         :: name
        
        ! Begin block
        call outFile % startBlock(self % getName())
        
        ! If collision clerk has map print map information
        if (allocated(self % map)) then
          call self % map % print(outFile)
        end if

        ! Write results.
        ! Get shape of result array
        if (allocated(self % map)) then
          resArrayShape = [1, self % map % binArrayShape()]
        else
          resArrayShape = [1]
        end if

    
        ! Start array
        name ='Res'
        call outFile % startArray(name, resArrayShape)

        ! Print results to the file
        do i = 1, product(resArrayShape)
          call mem % getResult(val, std, self % getMemAddress() + self % width * (i - 1) + 2)
          call outFile % addResult(val,std)

        end do

        call outFile % endArray()
        call outFile % endBlock()
      end subroutine print
    
    end module keffXClerk_class    



