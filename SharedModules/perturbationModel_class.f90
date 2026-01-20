module perturbationModel_class

  use numPrecision
  use genericProcedures,              only : fatalError, numToChar
  use dictionary_class,               only : dictionary
  implicit none


  type, public :: perturbationModel
    ! Define the attributes of the perturbation model here
    integer(shortInt)               :: perturbationType
    real(defReal), allocatable      :: densMagnitude(:)
    real(defReal)                   :: geometryMagnitude
    character(nameLen)              :: geometryType
    character(len=10), allocatable  :: perturbedIsotopes(:)
    integer(shortInt)               :: numPerturbations
  contains
    procedure :: init
    procedure :: kill
    procedure :: hasGeometricalPerturbation
  end type perturbationModel

  contains 

  function hasGeometricalPerturbation(self) result(hasGeomPert)
    class(perturbationModel), intent(in) :: self
    logical(defBool)                     :: hasGeomPert

    if (self % perturbationType == 1 .or. self % perturbationType == 3) then
        hasGeomPert = .true.
    else
        hasGeomPert = .false.
    end if
  end function hasGeometricalPerturbation

  subroutine init(self, dict)
    class(perturbationModel), intent(inout)       :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: keys
    class(dictionary), pointer                    :: tempDict
    integer(shortInt)                             :: i
    character(nameLen)                            :: tempChar  
    character(100), parameter                 :: Here = 'init (perturbationModel_class.f90)'

    ! Check if range of uncertain parameters is present
    if (dict % isPresent('type')) then
        call dict % get(tempChar, 'type')
        select case (trim(tempChar))
            case ('geometry')
                self % perturbationType = 1
            case ('density')
                self % perturbationType = 2
            case ('mixed')
                self % perturbationType = 3
            case default
                call fatalError(Here,  'Unknown perturbation type')
        end select
    else
        call fatalError(Here, 'Perturbation type not specified')
    end if
    self % numPerturbations = 0

    if (.not. self % perturbationType == 1) then
        ! Density perturbations
        if (.not. dict % isPresent('densityPerturbations')) then
            call fatalError(Here, 'No density perturbations provided in the dictionary')
        end if
        tempDict => dict % getDictPtr('densityPerturbations')
        call tempDict % keys(keys)

        allocate(self % densMagnitude(size(keys)))
        allocate(self % perturbedIsotopes(size(keys)))
        self % numPerturbations = size(keys)

        do i = 1, size(keys)
            self % perturbedIsotopes(i) = trim(keys(i))
            call tempDict % get(self % densMagnitude(i), keys(i))
        end do
    end if

    if (.not. self % perturbationType == 2) then
        ! Geometrical perturbations
        if (.not. dict % isPresent('geometryPerturbation')) then
            call fatalError(Here, 'No geometrical perturbation provided in the dictionary')
        end if
        tempDict => dict % getDictPtr('geometryPerturbation')

        ! Check if magnitude is provided
        if (tempDict % isPresent('magnitude')) then
            call tempDict % get(self % geometryMagnitude, 'magnitude')
        else
            call fatalError(Here, 'No magnitude provided for geometrical perturbation')
        end if

        if (tempDict % isPresent('type')) then
            call tempDict % get(self % geometryType, 'type')
            self % geometryType = trim(self % geometryType)
            if (self % geometryType /= 'isotropic' .and. &
                self % geometryType /= 'radial' .and. &
                self % geometryType /= 'axial') then
                call fatalError(Here, 'Unknown geometrical perturbation type: '//trim(self % geometryType))
            end if
        else
            call fatalError(Here, 'No type provided for geometrical perturbation')
        end if
        self % numPerturbations = self % numPerturbations + 1
    end if
  end subroutine init

  subroutine kill(self)
    class(perturbationModel), intent(inout) :: self

    if (allocated(self % densMagnitude)) deallocate(self % densMagnitude)
    if (allocated(self % perturbedIsotopes)) deallocate(self % perturbedIsotopes)
  end subroutine kill

end module perturbationModel_class