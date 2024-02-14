module fissionSourceClerk_class

   use numPrecision
   use tallyCodes
   use genericProcedures,          only : fatalError
   use dictionary_class,           only : dictionary
   use particle_class,             only : particle, particleState
   use particleDungeon_class,      only : particleDungeon
   use outputFile_class,           only : outputFile

   ! Basic tally modules
   use scoreMemory_class,          only : scoreMemory
   use tallyClerk_inter,           only : tallyClerk

   ! Tally Maps
   use tallyMap_inter,             only : tallyMap
   use tallyMapFactory_func,       only : new_tallyMap

   implicit none
   private

   !!
   !! Fission source estimator
   !! This is a prototype implementation of an analog estimator for the fission source
   !! Takes cycle start reports to generate cycle-wise fission source (after eventual population control has been applied)
   !! Contains only a single map for discretisation
   !! Scores fission source for a given number of cycles
   !!
   !! Notes:
   !!    ->
   !! Sample dictionary input:
   !!
   !!  clerkName {
   !!      type fissionSourceClerk;
   !!      map { <TallyMapDef> };
   !!  }
   !!
   type, public, extends(tallyClerk):: fissionSourceClerk
      private
      !! Map defining the discretisation
      class(tallyMap), allocatable                   :: map

   contains
      ! Procedures used during build
      procedure  :: init
      procedure  :: validReports
      procedure  :: getSize

      ! File reports and check status -> run-time procedures
      procedure  :: reportCycleStart

      ! Output procedures
      procedure  :: display
      procedure  :: print

      ! Deconstructor
      procedure  :: kill
   end type fissionSourceClerk

contains

   !!
   !! Initialise clerk from dictionary and name
   !!
   subroutine init(self, dict, name)
      class(fissionSourceClerk), intent(inout):: self
      class(dictionary), intent(in)             :: dict
      character(nameLen), intent(in)            :: name

      ! Assign name
      call self % setName(name)

      ! Read map
      call new_tallyMap(self % map, dict % getDictPtr('map'))

   end subroutine init

   !!
   !! Returns array of codes that represent diffrent reports
   !!
   function validReports(self) result(validCodes)
      class(fissionSourceClerk), intent(in)      :: self
      integer(shortInt), dimension(:), allocatable:: validCodes

      validCodes = [cycleStart_Code]

   end function validReports

   !!
   !! Return memory size of the clerk
   !!
   elemental function getSize(self) result(S)
      class(fissionSourceClerk), intent(in):: self
      integer(shortInt)                      :: S

      if (allocated(self % map)) S = self % map % bins(0)

   end function getSize


   !!
   !! Process cycle end
   !!
   subroutine reportCycleStart(self, start, mem)
      class(fissionSourceClerk), intent(inout):: self
      class(particleDungeon), intent(in)        :: start
      type(scoreMemory), intent(inout)          :: mem
      integer(shortInt)                         :: i, idx
      integer(longInt) :: addr
      type(particleState) :: state
      character(100), parameter :: Here = ' reportCycleStart (fissionSourceClerk_class.f90)'

      ! Loop over fission sites
      do i = 1, start % popSize()
         ! Get current particle
         state = start % get(i)

         ! Find bin index
         if (allocated(self % map )) then
            idx = self % map % map(state)
         else
            call fatalError(Here, 'Invalid mapping for fission source clerk')
         end if
         ! Calculate bin index address
         addr = self % getMemAddress() + idx - 1
         ! Score fission source contribution of current particle
         call mem % score(state % wgt, addr)
      end do

   end subroutine reportCycleStart

   !!
   !! Display convergance progress on the console
   !!
   subroutine display(self, mem)
      class(fissionSourceClerk), intent(in):: self
      type(scoreMemory), intent(in)    :: mem

      print *, 'fissionSourceClerk does not support display yet'

   end subroutine display

   !!
   !! Write contents of the clerk to output file
   !!
   subroutine print(self, outFile, mem)
      class(fissionSourceClerk), intent(in):: self
      class(outputFile), intent(inout):: outFile
      type(scoreMemory), intent(in)    :: mem
      integer(shortInt)                :: i
      integer(shortInt),dimension(:),allocatable :: resArrayShape
      character(nameLen)               :: name
      real(defReal) :: val, std

      ! Begin block
      call outFile % startBlock(self % getName())

      ! Print map information
      if (allocated(self % map)) then
         call self % map % print(outFile)
      end if

      ! Print fissionSource
      name = 'fissionSource'

      ! Get shape of result array
      if (allocated(self % map)) then
         resArrayShape = [1, self % map % binArrayShape()]
      end if

      ! Start array
      call outFile % startArray(name, resArrayShape)

      ! Write fission source results
      do i = 1, product(resArrayShape)
         call mem % getResult(val, std, self % getMemAddress() - 1 + i)
         call outFile % addResult(val, std)
      end do

      call outFile % endArray()
      call outFile % endBlock()

   end subroutine print

   !!
   !! Returns to uninitialised state
   !!
   elemental subroutine kill(self)
      class(fissionSourceClerk), intent(inout):: self

      if(allocated(self % map)) deallocate(self % map)

   end subroutine kill

end module fissionSourceClerk_class
