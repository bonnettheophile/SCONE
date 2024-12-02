module uncertainMap_class

  use numPrecision
  use universalVariables
  use genericProcedures,   only : fatalError
  use dictionary_class,    only : dictionary
  use grid_class,          only : grid
  use particle_class,      only : particleState
  use outputFile_class,    only : outputFile
  use tallyMap1D_inter,    only : tallyMap1D, kill_super => kill

  implicit none
  private

  !!
  !! Constructor
  !!
  interface uncertainMap
    module procedure uncertainMap_fromDict
  end interface

  !!
  !! Map that divides uncertain parameters into number of discrete bins
  !! Returns index 0 for elements outside division
  !!
  !! Private Members:
  !!   binBounds -> grid with bin boundaries
  !!   N -> Integer number of bins in the map
  !!
  !! Interface:
  !!   tallyMap interface
  !!   build -> build instance of spaceMap without dictionary
  !!
  !! NOTE: Behaviour of points exactly at the boundary between two bins is undefined.
  !!       They can be in either of two bins.
  !!
  !! Sample Dictionary Input:
  !!   structMap {
  !!     type uncertainMap;
  !!     grid lin;
  !!     min 0.0;
  !!     max 2.0;
  !!     N 20;
  !!   }
  !!
  !!   unstructMap {
  !!     type uncertainMap;
  !!     grid ustruct;
  !!     bins (0.0 0.2 0.5 0.8 1.0);
  !!   }
  !!
  type, public,extends(tallyMap1D) :: uncertainMap
    private
    type(grid)        :: binBounds
    integer(shortInt) :: N = 0
    integer(shortInt) :: dir

  contains
    ! Superclass interface implementaction
    procedure  :: init
    procedure  :: bins
    procedure  :: map
    procedure  :: getAxisName
    procedure  :: print
    procedure  :: kill

    ! Class specific procedures
    generic            :: build => build_fromGrid, build_structured
    procedure,private  :: build_fromGrid
    procedure,private  :: build_structured
  end type uncertainMap

contains

  !!
  !! Build from explicit grid of bin boundaries
  !!
  !! Args:
  !!   grid [in] -> Array of sorted ascending defReal values
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_fromGrid(self, grid)
    class(uncertainMap), intent(inout)         :: self
    real(defReal), dimension(:), intent(in) :: grid

    self % N = size(grid)-1
    call self % binBounds % init(grid)

  end subroutine build_fromGrid

  !!
  !! Build from min and max value, number of bins and direction
  !!
  !! Args:
  !!   mini [in] -> minumum value on the grid
  !!   maxi [in] -> maximum value on the grid
  !!   N [in] -> Number of bins in the grid
  !!   type [in] -> nameLen type of interpolation 'lin' or 'log'
  !!
  !! Errors:
  !!   None from here. Grid type is responsible for checking input consistency
  !!
  subroutine build_structured(self, mini, maxi, N, type)
    class(uncertainMap), intent(inout) :: self
    real(defReal), intent(in)       :: mini
    real(defReal), intent(in)       :: maxi
    integer(shortInt),intent(in)    :: N
    character(nameLen), intent(in)  :: type

    self % N = N
    call self % binBounds % init(mini, maxi, N, type)

  end subroutine build_structured

  !!
  !! Initialise from dictionary
  !!
  !! See tallyMap for specification
  !!
  subroutine init(self, dict)
    class(uncertainMap), intent(inout)        :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: str, type
    real(defReal)                          :: mini, maxi
    real(defReal),dimension(:),allocatable :: bins
    integer(shortInt)                      :: N
    character(100), parameter     :: Here = 'init (uncertainMap_class.f90)'

    if(.not.dict % isPresent('grid')) call fatalError(Here,"Keyword 'grid' must be present")

    ! Read grid definition keyword
    call dict % get(str,'grid')

    ! Choose approperiate definition
    select case(str)
      case('lin')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'lin'

        ! Initialise
        call self % build(mini, maxi, N, type)

      case('log')
        ! Read settings
        call dict % get(mini,'min')
        call dict % get(maxi,'max')
        call dict % get(N,'N')
        type = 'log'

        ! Initialise
        call self % build(mini, maxi, N, type)

      case('unstruct')
        ! Read settings
        call dict % get(bins,'bins')

        ! Initialise
        call self % build(bins)

      case default
        call fatalError(Here,"'grid' keyword must be: lin, log or usntruct")
    end select

    call dict % get(str, 'axis')

    select case(str)
      case('x')
        self % dir = 1
      case('y')
        self % dir = 2
      case('z')
        self % dir = 3
      case default
        call fatalError(Here, "'axis' keyword must be: x, y, or z")
    end select

  end subroutine init

  !!
  !! Return total number of bins in this division along dimension D
  !! For D=0 return all bins
  !!
  !! See tallyMap for specification
  !!
  elemental function bins(self, D) result(N)
    class(uncertainMap), intent(in)    :: self
    integer(shortInt), intent(in)   :: D
    integer(shortInt)               :: N

    if (D == 1 .or. D == 0) then
      N = self % N
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin.
  !!
  !! See tallyMap for specification
  !!
  elemental function map(self,state) result(idx)
    class(uncertainMap), intent(in)     :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    ! Find position on the grid
    idx = self % binBounds % search(state % X(self % dir))
    if (idx == valueOutsideArray) idx = 0

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(uncertainMap), intent(in) :: self
    character(nameLen)              :: name

    select case(self % dir)
      case(1)
        name = 'x'
      case(2)
        name = 'y'
      case(3)
        name = 'z'
    end select
  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification
  !!
  subroutine print(self,out)
    class(uncertainMap), intent(in)     :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Name the array
    name = trim(self % getAxisName()) //'Bounds'

    call out % startArray(name,[self % N,2])
    do i=1,self % N
      ! Print lower bin boundary
      call out % addValue(self % binBounds % bin(i))
    end do

    do i=1,self % N
      ! Print upper bin boundar
      call out % addValue(self % binBounds % bin(i+1))
    end do

    call out % endArray()

  end subroutine print

  !!
  !! Return instance of uncertainMap from dictionary
  !!
  !! Args:
  !!   dict[in] -> input dictionary for the map
  !!
  !! Result:
  !!   Initialised uncertainMap instance
  !!
  !! Errors:
  !!   See init procedure.
  !!
  function uncertainMap_fromDict(dict) result(new)
    class(dictionary), intent(in)          :: dict
    type(uncertainMap)                        :: new

    call new % init(dict)

  end function uncertainMap_fromDict

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uncertainMap), intent(inout) :: self

    call kill_super(self)

    ! Kill local
    call self % binBounds % kill()
    self % N = 0

  end subroutine kill

end module uncertainMap_class
