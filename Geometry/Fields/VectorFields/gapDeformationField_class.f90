module gapDeformationField_class

  use numPrecision
  use universalVariables
  use genericProcedures,        only : fatalError, numToChar, swap
  use field_inter,              only : field
  use deformationField_inter,   only : deformationField
  use vectorField_inter,        only : vectorField
  use coord_class,              only : coordList
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use box_class,                only : box
  use materialMenu_mod,         only : mm_matIdx => matIdx

  implicit none
  private

  !! Helper type for variable amplitude option
  type, public :: quadraticFunction
    real(defReal) :: a
    real(defReal) :: b
    real(defReal) :: c
    real(defReal) :: zLimit
    contains
      procedure :: init => init_quadraticFunction
      procedure :: eval
  end type quadraticFunction
  
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
  !!   variableAmplitude true; # If true, deformation amplitude is symmetric and quadratic wrt z. Max amplitude is rS,
  !!                           #  and zero deformation is applied at extremities of lattice in z. Default false.
  !!   # default 8.0; #
  !!
  !! }
  !!
  type, public, extends(deformationField) :: gapDeformationField
    type(quadraticFunction)    :: deformationAmplitude
    logical(defBool)           :: variableAmplitude = .false.
    real(defReal)              :: rI, rO, rMinus, rC
    real(defReal)              :: gapClosing  

  contains
    ! Superclass procedures
    procedure :: forward
    procedure :: backward
    procedure :: init 

  end type gapDeformationField

contains

  !! First implement helper type procedures
  !! Initialise quadratic function coefficients
  subroutine init_quadraticFunction(self, coeffs, zLimit)
    class(quadraticFunction), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: coeffs
    real(defReal), intent(in)               :: zLimit
    character(100), parameter :: Here = 'quadraticFunction % init (gapDeformationField_class.f90)'
    
    self % a = coeffs(1)
    self % b = coeffs(2)
    self % c = coeffs(3)
    self % zLimit = zLimit
  end subroutine init_quadraticFunction

  !! Evaluate quadratic function at given z value
  function eval(self, z, X) result(val)
    class(quadraticFunction), intent(in) :: self
    real(defReal), intent(in)           :: z
    real(defReal), intent(in), optional :: X
    real(defReal) :: val, X_value

    if (present(X)) then
      X_value = X
    else
      X_value = ONE
    end if

    if (z > X_value * self % zLimit) then
      val = ONE
    else if (z < ZERO) then
      val = ZERO
    else
      val = self % a * z**2 / X_value ** 2 + self % b * z / X_value + self % c
    end if

  end function eval

  !! Initialisation
  !!
  subroutine init(self, dict)
    class(gapDeformationField), intent(inout)          :: self
    class(dictionary), intent(in)                 :: dict
    type(dictionary)                              :: tempDict
    integer(shortInt)                             :: N, j, k, idx0
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

    ! Load deformation parameters
    call dict % getOrDefault(self % rI, 'innerAnnulus', 0.09_defReal)
    call dict % get(self % rO, 'outerAnnulus')
    call dict % get(self % rC, 'rCladding')
    call dict % get(self % rMinus, 'rMinus')
    call dict % get(self % variableAmplitude, 'variableAmplitude')
    call dict % getOrDefault(self % gapClosing, 'gapClosing', ZERO)

    if (self % variableAmplitude) then
      origin(1) = - ONE / (self % gapClosing ** 2)
      origin(2) = 2 / self % gapClosing
      origin(3) = ZERO
      call self % deformationAmplitude % init(origin, self % gapClosing)
      origin = ZERO
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

  !!
  !! Get value of the field at the co-ordinate point
  !!
  !!
  function forward(self, coords, X) result(r)
    class(gapDeformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(:), intent(in), optional :: X
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    real(defReal)                     :: zShift
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    ! In this case, we return the original coordinates (i.e. no deformation)
    if (localID == 0) then
      r = coords % lvl(1) % r
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
    zShift = self % pitch(3) / TWO

    delta = self % delta(localID)

    if (self % variableAmplitude .and. present(X)) then
      delta = delta * self % deformationAmplitude % eval(r_bar(3) + zShift, (X(1) + ONE) / TWO)
    else if (self % variableAmplitude) then
      delta = delta * self % deformationAmplitude % eval(r_bar(3) + zShift)
    end if

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    if (delta /= 0.0) then 
      if (r(1) > self % rMinus .and. r(1) < self % rO) then 
        r(1) = r(1) + delta / (self % rO - self % rMinus) * (r(1) - self % rMinus)
      elseif (r(1) > self % rO .and. r(1) < self % rC) then 
        r(1) = r(1) + delta / (self % rC - self % rO) * (self % rC - r(1))
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
  function backward(self, coords, X) result(r)
    class(gapDeformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(:), intent(in), optional :: X
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    real(defReal)                     :: zShift
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    ! In this case, we return the original coordinates (i.e. no deformation)
    if (localID == 0) then
      r = coords % lvl(1) % r
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
    zShift = self % pitch(3) / TWO

    delta = self % delta(localID)

    if (self % variableAmplitude .and. present(X)) then
      delta = delta * self % deformationAmplitude % eval(r_bar(3) + zShift, (X(1) + ONE) / TWO)
    else if (self % variableAmplitude) then
      delta = delta * self % deformationAmplitude % eval(r_bar(3) + zShift)
    end if

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    if (delta /= 0.0) then 
      if (r(1) > self % rMinus .and. r(1) < self % rO + delta) then 
        r(1) = (r(1) + delta / (self % rO - self % rMinus) * self % rMinus) / (1 + delta / (self % rO - self % rMinus))
      elseif (r(1) > self % rO + delta .and. r(1) < self % rC) then 
        r(1) = (r(1) - delta * self % rC / (self % rC - self % rO)) / (1 - delta / (self % rC - self % rO))
      end if
    end if
        
    r = self % xTransform(r)

    ! Revert to original coordinates
    r = r + self % corner + (ijk - HALF) * self % pitch

  end function backward
  

end module gapDeformationField_class
