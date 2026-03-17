module deformationField_class

  use numPrecision
  use universalVariables
  use genericProcedures,        only : fatalError, numToChar, swap
  use field_inter,              only : field
  use vectorField_inter,        only : vectorField
  use coord_class,              only : coordList
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use box_class,                only : box
  use materialMenu_mod,         only : mm_matIdx => matIdx

  implicit none
  private
  
  !!
  !! Public Pointer Cast
  !!
  public :: deformationField_TptrCast

  integer(shortInt), parameter :: ALL_MATS = -1

  !!
  !! Piecewise constant field constructed from a lattice-like grid. 
  !! Values of the field are piecewise constant.
  !!
  !! Similar to a Cartesian lattice. Centre is placed at origin.
  !! delta is a deformation amplitude array
  !!
  !! Example dictionary:
  !!
  !! myField {
  !!   type deformationField;
  !!   origin (x0 y0 z0);
  !!   shape (Nx Ny Nz);
  !!   pitch (Px Py Pz);
  !!   ! Make up to size Nx * Ny * Nz, ascending first in x, then y, then z
  !!   delta (
  !!    100 92 3.14 ...
  !!   ); 
  !!   rS 0.54; # Radius of fuel pin to be deformed
  !!   r0 0.63; # Outer radius of deformation region. Should be larger than rS.
  !!   # default 8.0; #
  !!
  !! }
  !!
  type, public, extends(vectorField) :: deformationField
    private
    real(defReal), dimension(3)               :: pitch = ZERO
    integer(shortInt), dimension(3)           :: sizeN = 0
    integer(shortInt)                         :: N = 0
    real(defReal), dimension(3)               :: corner = ZERO
    real(defReal), dimension(3)               :: a_bar  = ZERO
    real(defReal), dimension(:), allocatable  :: delta
    real(defReal)                             :: rS 
    real(defReal)                             :: r0
    type(box)                                 :: outline
    
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: at
    procedure :: atP
    procedure :: forward
    procedure :: backward
    procedure :: rTransform
    procedure :: xTransform
    procedure :: kill

    procedure, private :: getLocalID
  end type deformationField

contains

  !!
  !! Initialisation
  !!
  subroutine init(self, dict)
    class(deformationField), intent(inout)          :: self
    class(dictionary), intent(in)                 :: dict
    type(dictionary)                              :: tempDict
    integer(shortInt)                             :: N, i, j, k, idx0
    integer(shortInt), dimension(:), allocatable  :: tempI
    real(defReal), dimension(:), allocatable      :: temp
    real(defReal), dimension(3)                   :: origin
    real(defReal), dimension(:,:), allocatable    :: tempMap
    character(100), parameter :: Here = 'init (deformationField_class.f90)'
    
    ! Load pitch
    call dict % get(temp, 'pitch')
    N = size(temp)

    if (N /= 3) then
      call fatalError(Here, 'Pitch must have size 3. Has: '//numToChar(N))
    end if
    self % pitch = temp
    
    ! Load origin
    call dict % get(temp, 'origin')
    N = size(temp)

    if (N /= 3) then
      call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(N))
    end if
    origin = temp

    ! Load Size
    call dict % get(tempI, 'shape')
    N = size(tempI)

    if (N /= 3) then
      call fatalError(Here, 'Shape must have size 3. Has: '//numToChar(N))
    else if (any(tempI < 0)) then
      call fatalError(Here, 'Shape contains -ve entries')
    end if
    self % sizeN = tempI

    ! Detect reduced Z dimension
    if (self % sizeN(3) == 0) then
      self % sizeN(3) = 1
      self % pitch(3) = TWO * INF
    end if

    ! Check X & Y for 0 size
    if (any( self % sizeN == 0)) call fatalError(Here, 'Shape in X and Y axis cannot be 0.')

    ! Check for invalid pitch
    if (any(self % pitch < 10 * SURF_TOL)) then
     call fatalError(Here, 'Pitch size must be larger than: '//numToChar( 10 * SURF_TOL))
   end if

    ! Calculate halfwidth and corner
    self % a_bar = self % pitch * HALF - SURF_TOL
    self % corner = origin -(self % sizeN * HALF * self % pitch)

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', origin)
    call tempDict % store('halfwidth', abs(self % corner - origin))
    call self % outline % init(tempDict)

    ! Size field value array
    self % N = product(self % sizeN)
    allocate(self % delta(self % N + 1))

    ! Read field values for each material
    idx0 = 0
    call dict % get(temp, 'delta')

      ! Flip array up-down for more natural input
      ! Reshape into rank 2 array
    tempMap = reshape(temp, [self % sizeN(1), self % sizeN(2) * self % sizeN(3)])
    N = size(tempMap, 2)
    do j = 1, N/2
      call swap(tempMap(:,j), tempMap(:,N - j + 1))
    end do

      ! Build fill array
      N = size(tempMap, 1)
      do j = 1, size(tempMap, 2)
        do k = 1, N
          self % delta(idx0 + k + (j-1) * N) = tempMap(k, j)
        end do
      end do

  end subroutine init

  function rTransform(self, x) result(r)
    class(deformationField), intent(in) :: self
    real(defReal), intent(in)              :: x(3)
    real(defReal)                          :: r(3)

    r(1) = sqrt(sum(x(1:2)**2))
    r(2) = atan2(x(2), x(1))
    r(3) = x(3)

  end function rTransform

  function xTransform(self, r) result(x)
    class(deformationField), intent(in) :: self
    real(defReal), intent(in)              :: r(3)
    real(defReal)                          :: x(3)

    x(1) = r(1)*cos(r(2))
    x(2) = r(1)*sin(r(2))
    x(3) = r(3)

  end function xTransform

  !!
  !! Get value of the field at the co-ordinate point
  !!
  !!
  function forward(self, coords) result(r)
    class(deformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta, d
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    if (localID == 0) then
      d = self % outline % distance(coords % lvl(1) % r, coords % lvl(1) % dir)
      return

    end if

    ! Compute ijk of localID
    temp = localID - 1

    base = temp / self % sizeN(1)
    ijk(1) = temp - self % sizeN(1) * base + 1

    temp = base
    base = temp / self % sizeN(2)
    ijk(2) = temp - self % sizeN(2) * base + 1

    ijk(3) = base + 1

    ! Find position wrt lattice cell centre
    ! Need to use localID to properly handle under and overshoots
    r_bar = coords % lvl(1) % r - self % corner
    r_bar = r_bar - (ijk - HALF) * self % pitch

    delta = self % delta(localID)

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    if (delta /= 0.0) then 
      if (r(1) < self % rS) then 
        r(1) = r(1) + delta / self % rS * r(1)
      elseif (r(1) > self % rS .and. r(1) < self % r0) then 
        r(1) = r(1) + delta / (self % rS - self % r0) * (r(1) - self % r0)
      end if
    end if
        
    r = self % xTransform(r)

    ! Revert to original coordinates
    r = r + self % corner + (ijk - HALF) * self % pitch

  end function forward

    !!
  !! Get value of the field at the co-ordinate point
  !!
  !!
  function backward(self, coords) result(r)
    class(deformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta, d
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    if (localID == 0) then
      d = self % outline % distance(coords % lvl(1) % r, coords % lvl(1) % dir)
      return

    end if

    ! Compute ijk of localID
    temp = localID - 1

    base = temp / self % sizeN(1)
    ijk(1) = temp - self % sizeN(1) * base + 1

    temp = base
    base = temp / self % sizeN(2)
    ijk(2) = temp - self % sizeN(2) * base + 1

    ijk(3) = base + 1

    ! Find position wrt lattice cell centre
    ! Need to use localID to properly handle under and overshoots
    r_bar = coords % lvl(1) % r - self % corner
    r_bar = r_bar - (ijk - HALF) * self % pitch

    delta = self % delta(localID)

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    if (delta /= 0.0) then 
      if (r(1) < self % rS + delta) then 
        r(1) = r(1) / (1 + delta / self % rS)
      elseif (r(1) > self % rS + delta .and. r(1) < self % r0) then 
        r(1) = (r(1) + delta * self % r0 / (self % rS - self % r0)) / (1 + delta / (self % rS - self % r0))
      end if
    end if
        
    r = self % xTransform(r)

    ! Revert to original coordinates
    r = r + self % corner + (ijk - HALF) * self % pitch

  end function backward
  
  !!
  !! Get value of the field at the particle's location
  !!
  !! See pieceConstantField for details
  !!
  function at(self, coords) result(val)
    class(deformationField), intent(in) :: self
    class(coordList), intent(in)        :: coords
    real(defReal), dimension(3)         :: val

    val = 0.0

  end function at

  !!
  !! Get value of the field at the particle's location
  !!
  !! See pieceConstantField for details
  !!
  function atP(self, p) result(val)
    class(deformationField), intent(in) :: self
    class(particle), intent(in)       :: p
    real(defReal), dimension(3)       :: val

    val = self % at(p % coords)

  end function atP
    
  !!
  !! Clean-up
  !!
  elemental subroutine kill(self)
    class(deformationField), intent(inout) :: self
    
    self % pitch = ZERO
    self % sizeN = 0
    self % corner = ZERO
    self % a_bar  = ZERO
    call self % outline % kill()

  end subroutine kill

  !!
  !! Find the local integer ID in the field given position and direction
  !!
  pure function getLocalID(self, r, u) result(localID)
    class(deformationField), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: localID
    real(defReal), dimension(3)             :: r_bar
    integer(shortInt), dimension(3)         :: ijk
    integer(shortInt)                       :: i, inc
          
    ijk = floor((r - self % corner) / self % pitch) + 1

    ! Get position wrt middle of the lattice cell
    r_bar = r - self % corner - ijk * self % pitch + HALF * self % pitch

    ! Check if position is within surface tolerance
    ! If it is, push it to next cell
    do i = 1, 3
      if (abs(r_bar(i)) > self % a_bar(i) .and. r_bar(i)*u(i) > ZERO) then

        ! Select increment. Ternary expression
        if (u(i) < ZERO) then
          inc = -1
        else
          inc = 1
        end if

        ijk(i) = ijk(i) + inc
      end if
    end do

    if (any(ijk <= 0 .or. ijk > self % sizeN)) then ! Point is outside lattice
      localID = 0

    else
      localID = ijk(1) + self % sizeN(1) * (ijk(2)-1 + self % sizeN(2) * (ijk(3)-1))

    end if


  end function getLocalID
    
  !!
  !! Cast field pointer to deformationField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of deformationField
  !!   Pointer to source if source is deformationField type
  !!
  pure function deformationField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(deformationField), pointer     :: ptr

    select type (source)
      type is (deformationField)
        ptr => source

      class default
        ptr => null()
    end select

  end function deformationField_TptrCast

end module deformationField_class
