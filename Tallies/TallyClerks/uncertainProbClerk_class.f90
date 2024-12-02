module uncertainProbClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use dictionary_class,      only : dictionary
    use genericProcedures,     only : fatalError
    use particle_class,        only : particle, particleState
    use particleDungeon_class, only : particleDungeon
    use outputFile_class,      only : outputFile
    use legendrePoly_func,     only : evaluateLegendre
    use polynomial_class,      only : polynomial

    use scoreMemory_class,     only : scoreMemory
    use tallyResult_class,     only : tallyResult, tallyResultEmpty, polyResult
    use tallyClerk_inter,      only : tallyClerk, kill_super => kill

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap

    use linearAlgebra_func
  
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
    !! myuncertainProbClerk {
    !!    type uncertainProbClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: uncertainProbClerk
      private
        class(tallyMap), allocatable             :: map
        integer(shortInt)                        :: fitOrder     ! Order of polyfit for pdf
        real(defReal), allocatable               :: histogram(:)
        real(defReal), allocatable               :: binCentre(:)
        real(defReal), allocatable               :: fitCoeff(:,:)
        logical(defBool)                         :: inactive = .false.


    contains
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize

      procedure :: reportCycleEnd

      procedure :: print
      procedure :: display
      procedure :: getResult
    end type uncertainProbClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(uncertainProbClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      integer(shortInt)                       :: i
      real(defReal)                           :: dx
      character(100),parameter :: Here = 'init (uncertainProbClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      call dict % getOrDefault(self % inactive, 'inactive', .false.)

      ! Order keyword must be present
      if (dict % isPresent('fitOrder')) then
        call dict % get(self % fitOrder, ' fitOrder')
        allocate(self % fitCoeff(1, self % fitOrder + 1))
        self % fitCoeff = ZERO
      else 
        call fatalError(Here, "fitOrder must by provided") 
      end if

      ! Map is for following gpc coefficients
      if( dict % isPresent('map')) then
        call new_tallyMap(self % map, dict % getDictPtr('map'))
        allocate(self % histogram(product(self % map % binArrayShape())))
        allocate(self % binCentre(product(self % map % binArrayShape())))

        self % histogram = ZERO

        dx = TWO / size(self % histogram)
        self % binCentre(1) = - ONE + dx / TWO

        ! Initialize binCentre, assume interval is [-1,1]
        do i = 2, size(self % binCentre)
          self % binCentre(i) = self % binCentre(i-1) + dx 
        end do
      else
        call fatalError(Here, 'No map has been defined for polyChaos')
      end if


    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(uncertainProbClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % fitCoeff)) deallocate(self % fitCoeff)
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      if (allocated(self % histogram)) deallocate(self % histogram)
      if (allocated(self % binCentre)) deallocate(self % binCentre)

      self % fitOrder = 0

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(uncertainProbClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [ cycleEnd_CODE ]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(uncertainProbClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      S = self % fitOrder + 1
    end function getSize

    !! Process end of cycle and compute chaotic coefficients 
    !! Assuming we have a single parameter
    !! 
    subroutine reportCycleEnd(self, end, mem)
      class(uncertainProbClerk), intent(inout)     :: self
      class(particleDungeon), intent(in)          :: end
      type(scoreMemory), intent(inout)            :: mem
      type(particleState)                         :: state
      type(particle)                              :: p
      integer(shortInt)                           :: i, binIdx
      real(defReal), dimension(size(self % binCentre))                      :: x, b
      real(defReal), dimension(size(self % binCentre), self % fitOrder+1)   :: A 
      character(100),parameter :: Here = 'reportCycleEnd (uncertainProbClerk.f90)'

      ! Reinialize histogram array
      self % histogram = ZERO
      do i = 1, end % popSize()
        p = end % get(i)
        state = p 

        ! Find bin index
        if (allocated(self % map)) then
          binIdx = self % map % map(state)
        else
          binIdx = 1
        end if
        ! Return if invalid bin index
        if (binIdx == 0) return
        
        ! Fill histogram of particles wrt their uncertain parameter
        self % histogram(binIdx) = self % histogram(binIdx) + state % wgt
      end do

      ! Set x array for linear fitting
      do i = 1, size(A, 2)
        A(:,i) = self % binCentre**(i-1)
      end do 

      ! Set y for linear fitting
      b = self % histogram
      ! Perform least square linear fitting using LAPACK 
      call solveLeastSquare(A, x, b)
      ! Save fit results
      if (.not. self % inactive) then
        self % fitCoeff(1,:) = self % fitCoeff(1,:) + (x(1:self % fitOrder+1) - self % fitCoeff(1,:)) / (mem % cycles+1)
      else
        self % fitCoeff(1,:) = x(1:self % fitOrder+1) 
      end if
      call kill_linearAlgebra()

    end subroutine reportCycleEnd

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(uncertainProbClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, "Fitting order : ", self % fitOrder, "Fitting coeffs : ", self % fitCoeff(1,:)
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(uncertainProbClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
    
        ! Do nothing
    
      end subroutine print
    
  !! 
  !! Get fitting model for end of generation probability distribution of X
  !!

    pure subroutine getResult(self, res, mem)
      class(uncertainProbClerk), intent(in)              :: self
      class(tallyResult), allocatable, intent(inout)    :: res
      type(scoreMemory), intent(in)                     :: mem

      allocate(res, source= polyResult(self % fitCoeff))

    end subroutine
      
  end module uncertainProbClerk_class    



