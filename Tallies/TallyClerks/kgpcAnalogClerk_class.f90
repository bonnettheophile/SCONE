module kgpcAnalogClerk_class
    use numPrecision
    use tallyCodes
    use endfConstants
    use genericProcedures,          only : fatalError, charCmp
    use dictionary_class,           only : dictionary
    use particle_class,             only : particle
    use particleDungeon_class,      only : particleDungeon
    use outputFile_class,           only : outputFile
    use legendrePoly_func
  
    ! Nuclear Data Interfaces
    use nuclearDataReg_mod,         only : ndReg_get => get
    use nuclearDatabase_inter,      only : nuclearDatabase
    use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast
    use neutronXSPackages_class,    only : neutronMacroXSs
  
    ! Tally Interfaces
    use scoreMemory_class,          only : scoreMemory
    use tallyResult_class,          only : tallyResult, tallyResultEmpty
    use tallyClerk_inter,           only : tallyClerk, kill_super => kill
    use keffAnalogClerk_class,      only : keffResult
  
    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap
  
    implicit none
    private
  
    type, public, extends(tallyClerk) :: kgpcAnalogClerk
      private
        ! Holder for gpc order
        integer(shortInt)           :: P
        ! Holder for number of uncertain parameters
        integer(shortInt)           :: Q
  
    contains
      ! Duplicate interface of the tallyClerk
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize
  
      ! File reports and check status -> run-time procedures
      procedure :: reportCycleEnd
  
      ! Output procedures
      procedure :: display
      procedure :: print
    end type kgpcAnalogClerk
      
  contains
  
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details 
  !!
  
  subroutine init(self, dict, name)
      class(kgpcAnalogClerk), intent(inout)     :: self
      class(dictionary), intent(in)               :: dict
      character(nameLen), intent(in)              :: name
      character(100), parameter                   :: Here = 'init (kgpcAnalogClerk_class.f90)'
  
      ! Set name
      call self % setName(name)
  
      ! Load gpc poynomial order
      if (dict % isPresent('order')) then
        call dict % get(self % Q, 'order')
      else 
        call fatalError(Here, 'Must include gpc polynomial order')
      end if
  
      ! Load number of uncertain parameters
      if (dict % isPresent('params')) then
        call dict % get(self % P, 'params')
      else 
        call fatalError(Here, 'Must include number of uncertain parameters')
      end if
  
  end subroutine init
  
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    !! See tallyClerk_inter for details
    !!
  function validReports(self) result(validCodes)
      class(kgpcAnalogClerk), intent(in)      :: self
      integer(shortInt), dimension(:), allocatable :: validCodes
  
      validCodes = [cycleEnd_CODE]
  
  end function validReports
  
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!
  elemental function getSize(self) result(S)
    class(kgpcAnalogClerk), intent(in) :: self
    integer(shortInt)                    :: S
    S = (self % Q + 1) ** self % P
  end function getSize

    subroutine reportCycleEnd(self, end, mem)
        class(kgpcAnalogClerk), intent(inout)                    :: self
        class(particleDungeon), intent(in)                       :: end
        type(scoreMemory), intent(inout)                         :: mem
        real(defReal)                                            :: k_est
        integer(shortInt), dimension(self % P + 1)               :: ind
        integer(shortInt)                                        :: max_ind
        integer(shortInt)                                        :: i,j,k
        real(defReal), dimension(self % P, self % Q + 1)         :: pol_eval
        real(defReal)                                            :: pol_product
        integer(longInt)                                         :: addr
    
        do j = 1, end % popSize()
            do i = 1, size(pol_eval,1)
                pol_eval(i,:) = evaluateLegendre(self % Q, end % prisoners(j) % oldX(1,i))
            end do
            
            ! Get history-wise k_eff
            k_est = end % prisoners(j) % k_eff

            ! Set max index to polynomial order + 1
            max_ind = self % Q + 1
            k = 1
            ind = ONE
        
            ! Dynamic depth loop for scoring all gpc coefficients
            do while (ind(size(ind)) == 1 )
                pol_product = ONE
                addr = ZERO
                do i = 1, size(pol_eval, 1)
                ! Compute product of legendre polynomials at given order subset and X
                pol_product = (2 * (ind(i)-1) + 1) * pol_product * pol_eval(i, ind(i))
                ! Compute adress for scoring given coefficient
                addr = addr + (ind(i) - 1) * (max_ind) ** (i - 1)
                end do
                ! Project keff onto gpc basis to get coefficients + wgt of MC integration
                pol_product = pol_product * k_est / end % popSize()
                ! Get adress in memory and score
                addr = addr + self % getMemAddress()
                call mem % score(pol_product, addr)
                ! Update indices of nested loop
                ind(1) = ind(1) + 1;
                do while(ind(k) == max_ind + 1)
                ind(k) = 1
                k = k + 1
                ind(k) = ind(k) + 1
                if (ind(k) /= max_ind + 1) k = 1
                end do
            end do
        end do
    
      end subroutine reportCycleEnd
  
    !!
    !! Display convergance progress on the console
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine display(self, mem)
      class(kgpcAnalogClerk), intent(in)  :: self
      type(scoreMemory), intent(in)         :: mem
  
      print *, 'collisionClerk does not support display yet'
  
    end subroutine display
  
    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(kgpcAnalogClerk), intent(inout) :: self
  
      ! Superclass
      call kill_super(self)
      
      self % P = 0
      self % Q = 0
    end subroutine kill
  
    !!
    !! Write contents of the clerk to output file
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine print(self, outFile, mem)
      class(kgpcAnalogClerk), intent(in)       :: self
      class(outputFile), intent(inout)           :: outFile
      type(scoreMemory), intent(in)              :: mem
      real(defReal)                              :: val, std
      integer(shortInt)                          :: i
      integer(shortInt),dimension(self % P)      :: resArrayShape
      character(nameLen)                         :: name
     
      ! Begin block
      call outFile % startBlock(self % getName())
  
      ! Write results.
      ! Get shape of result array
      resArrayShape = self % Q + 1
  
      ! Start array
      name ='Res'
      call outFile % startArray(name, resArrayShape)
  
      ! Print results to the file
      do i=1,product(resArrayShape)
        call mem % getResult(val, std, self % getMemAddress() - 1 + i)
        call outFile % addResult(val,std)
  
      end do
  
      call outFile % endArray()
      call outFile % endBlock()
  
    end subroutine print
  end module kgpcAnalogClerk_class