module polynomial_class

    use numPrecision
    use genericProcedures,      only : fatalError

    implicit none
    private

    type, public :: polynomial
        private
        real(defReal), allocatable :: coefficients(:,:)
        integer(shortInt)          :: order = -1
        integer(shortInt)          :: dimension = -1 
    contains

        procedure :: buildND
        procedure :: build1D
        generic   :: build => build1D, buildND
        generic   :: evaluate => evaluateVector, evaluateVectorArray
        procedure :: deg
        procedure :: d
        procedure :: evaluateVector
        procedure :: evaluateVectorArray
        procedure :: kill
    end type polynomial

contains

    ! Initializes the polynomial with given coefficients for ND polynomials
    subroutine buildND(self, coeffs)
        class(polynomial), intent(inout) :: self
        real(defReal), intent(in)        :: coeffs(:,:)

        if (allocated(self % coefficients)) deallocate(self % coefficients)

        allocate(self % coefficients, source=coeffs)
        self % dimension = size(self % coefficients, dim= 1)
        self % order = size(self % coefficients, dim= 2) - 1
    end subroutine buildND

    ! Initializes the polynomial with given coefficients for 1D polynomials
    subroutine build1D(self, coeffs)
        class(polynomial), intent(inout) :: self
        real(defReal), intent(in)        :: coeffs(:)

        if (allocated(self % coefficients)) deallocate(self % coefficients)

        allocate(self % coefficients(1, size(coeffs)))
        self % coefficients(1,:) = coeffs
        self % dimension = size(self % coefficients, dim= 1)
        self % order = size(self % coefficients, dim= 2) - 1
    end subroutine build1D

    ! Evaluates ND polynomial at a given x-vector
    ! Args:
    !   x -> 1d array with size of polynomial dimension
    ! Return a single value
    function evaluateVector(self, x) result(val)
        class(polynomial), intent(in) :: self
        real(defReal), intent(in)     :: x(:)
        real(defReal)                 :: val
        integer(shortInt)             :: i, j
        character(100),parameter      :: Here = 'evaluateVector (polynomial_class.f90)'

        if (size(x) /= self % dimension) call fatalError(Here, "Invalid dimension of x")

        val = ZERO
        do i = 1, self % dimension
            do j = 1, self % order + 1
                val = val + self % coefficients(i,j) * x(i)**(j - 1)
            end do
        end do
    end function evaluateVector

    ! Evaluates ND polynomial at for an array of x-vectors
    ! Args:
    !   x -> 2d array with dim 1 = polynomial dimension
    !                      dim 2 = number of points to evaluate
    ! Return a 1D array
    function evaluateVectorArray(self, x) result(val)
        class(polynomial), intent(in) :: self
        real(defReal), intent(in)     :: x(:,:)
        real(defReal)                 :: val(size(x, dim=2))
        integer(shortInt)             :: i, j
        character(100),parameter :: Here = 'evaluateVectorArray (polynomial_class.f90)'

        if (size(x, dim=1) /= self % dimension) call fatalError(Here, "Invalid dimension of x")

        

        val = ZERO
        do i = 1, self % dimension
            do j = 1, self % order + 1
                val = val + self % coefficients(i,j) * x(i,:)**(j - 1)
            end do
        end do
    end function evaluateVectorArray

    ! Returns the polynomial order
    pure function deg(self) result(val)
        class(polynomial), intent(in)   :: self
        integer(shortInt)               :: val

        val = self % order
    end function deg 

    pure function d(self) result(val)
        class(polynomial), intent(in)   :: self
        integer(shortInt)               :: val

        val = self % dimension
    end function d

    ! Deallocates the polynomial and resets it
    elemental subroutine kill(self)
        class(polynomial), intent(inout) :: self

        if (allocated(self % coefficients)) deallocate(self % coefficients)
        self % order = -1
        self % dimension = -1
    end subroutine kill

end module polynomial_class