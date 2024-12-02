module implicitChaosClerk_class

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
    use genericProcedures,     only : binarySearch, quickSort
  
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
    !!    type implicitChaosClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: implicitChaosClerk
      private
        integer(shortInt)                        :: P            ! Order of gpc model
        real(defReal)                            :: startPop
      
    contains
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize

      procedure :: reportCycleStart
      procedure :: reportHist

      procedure :: print
      procedure :: display

    end type implicitChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(implicitChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      character(100),parameter :: Here = 'init (implicitChaosClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      ! Order keyword must be present
      if (dict % isPresent('order')) then
        call dict % get(self % P, ' order')
      else 
        call fatalError(Here, "Order must by provided") 
      end if

      self % startPop = ZERO
    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(implicitChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(implicitChaosClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [hist_CODE, cycleStart_CODE]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(implicitChaosClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      S = self % P + 1
    end function getSize

    subroutine reportCycleStart(self, start, mem)
      class(implicitChaosClerk), intent(inout) :: self
      class(particleDungeon), intent(in)      :: start
      type(scoreMemory), intent(inout)        :: mem

      self % startPop = start % popWeight()

    end subroutine reportCycleStart
  
  
    !!
    !! Process history report
    !! Gets fate code from the particle
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportHist(self, p, xsData, mem)
      class(implicitChaosClerk), intent(inout) :: self
      class(particle), intent(in)             :: p
      class(nuclearDatabase),intent(inout)    :: xsData
      type(scoreMemory), intent(inout)        :: mem
      real(defReal)                           :: score, val
      real(defReal), dimension(self % P + 1)  :: legendrePol
      integer(shortInt)                       :: j


      if (p % fate /= leak_FATE) then
        val = p % X(p % gpcPert)
        ! Evaluate Legendre polynomials up to right order
        legendrePol = evaluateLegendre(self % P, val) 
        do j = 1, self % P + 1
          score = (2*(j-1) + 1) * legendrePol(j) * p % w / self % startPop
          call mem % score(score, self % getMemAddress() + j - 1)
        end do
      end if
    end subroutine reportHist

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(implicitChaosClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'implicitChaosClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(implicitChaosClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: i
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        character(nameLen)                         :: name
    
        ! Begin block
        call outFile % startBlock(self % getName())
    
        ! Write results.
        ! Get shape of result array
        resArrayShape = [self % getSize()]
    
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
    
    end module implicitChaosClerk_class    



