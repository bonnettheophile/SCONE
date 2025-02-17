module collisionProcessorFactory_func

   use numPrecision
   use genericProcedures, only : fatalError
   use dictionary_class,  only : dictionary

   ! Abstract interface
   use collisionProcessor_inter, only : collisionProcessor

   ! Implementation
   use neutronCEstd_class, only : neutronCEstd
   use neutronCEimp_class, only : neutronCEimp
   use neutronCEblm_class, only : neutronCEblm
   use neutronCEbls_class, only : neutronCEbls
   use neutronMGstd_class, only : neutronMGstd
   use neutronMGbls_class, only : neutronMGbls

   implicit none
   private

   public:: new_collisionProcessor

   ! List that contains all accaptable types of collisionProcessors
   ! It is printed if type was unrecognised
   ! NOTE:
   ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
   character(nameLen), dimension(*), parameter:: AVALIBLE_collisionProcessors = [ 'neutronCEstd',&
      'neutronCEimp',&
      'neutronMGstd',&
      'neutronMGbls',&
      'neutronCEbls',&
      'neutronCEblm']

contains

   !!
   !! Allocate new allocatable collisionProcessor to a specific type
   !! If new is allocated it deallocates it
   !!
   subroutine new_collisionProcessor(new, dict)
      class(collisionProcessor), allocatable, intent(inout):: new
      class(dictionary), intent(in)                        :: dict
      character(nameLen)                                   :: type
      character(100), parameter      :: Here = 'new_collisionProcessor (collisionProcessorFactory_func.f90)'

      ! Deallocate new if allocated
      if(allocated(new)) deallocate(new)

      ! Obtain string that specifies type to be built
      call dict % get(type, 'type')

      ! Allocate approperiate subclass of collisionProcessor
      select case(type)
       case('neutronCEstd')
         allocate(neutronCEstd:: new)

       case('neutronCEimp')
         allocate(neutronCEimp:: new)

       case('neutronCEbls')
         allocate(neutronCEbls:: new)

       case('neutronCEblm')
         allocate(neutronCEblm:: new)

       case('neutronMGstd')
         allocate(neutronMGstd:: new)

       case('neutronMGbls')
         allocate(neutronMGbls:: new)

       case default
         print *, AVALIBLE_collisionProcessors
         call fatalError(Here, 'Unrecognised type of collisionProcessor: ' // trim(type))

      end select

      ! Initialise new processor
      call new % init(dict)

   end subroutine new_collisionProcessor

end module collisionProcessorFactory_func
