module legendrePoly_func
  !! This module contains functions related to evaluating and sampling Legendre Polynomials
  !! All of the methods here are based on:
  !!   "Monte Carlo Particle Transport Methods: Neutron and Photon Calculations"
  !!    by I. Lux and L. Koblinger
  !!    ISBN 978-1-315-89573
  !!

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG

  implicit none
  private

  ! Public Interface
  public :: sampleLegendre
  public :: evaluateLegendre

  ! Generic Interfaces

  !! Evaluate value of Legendre polynomial of order n at given point
  !! Might include different methods to compute the n-th Legendre polynomial
  interface evaluateLegendre
    module procedure evaluateLegendre_Pn
  end interface evaluateLegendre

  !!
  !! Samples diffrent Legendre PDFs given number of coefficients and RNG
  !!   Diffrent Orders can use diffrent methods and have some limits on supported coefficients
  !!   Consult particular implementations for details
  !!

  interface sampleLegendre
    module procedure sampleLegendre_P1
  end interface sampleLegendre

contains
  
  !! 
  !! Function to compute recursively the n-th Legendre polynomial
  !!
  !! Return the vector for Legendre polynomials evaluated at "value" up to order n
  !!
  !! Args: 
  !!   n [in] -> largest order of sought Legendre polynomials
  !!   value [in] -> point at which to evaluate Legendre polynomials
  !!
  !! Return:
  !!   array with the value of the polynomials evaluated at value
  !! 
  !! Error:
  !!   fatalError if abs(value) > 1
  function evaluateLegendre_Pn(n, value) result(val_vec)
    integer(shortInt), intent(in)         :: n
    real(defReal), intent(in)             :: value
    real(defReal), dimension(max(2,n+1))    :: val_vec
    integer(shortInt)                     :: i
    character(100), parameter :: Here = 'evaluateLegendre_Pn ( legendrePoly_func.f90)'

    ! Checks if value is valid
    if (abs(value) > 1) call fatalError(Here, 'Legendre polynomials are defined in [-1,1]')

    ! Return at least the first two polynomials, as there is no point in having n = 0
    val_vec(1) = ONE
    val_vec(2) = value

    ! Computes up to n-th using Bonnet recurrence formula
    do i = 2, n
      val_vec(i+1) = ((2*i-1)*value*val_vec(i)-(i-1)*val_vec(i-1))/(i)
    end do
    
  end function evaluateLegendre_Pn

  

  !!
  !! Simple function to sample P1 Legendre PDF. Uses methods from Lux & Koblinger
  !!
  function sampleLegendre_P1(P1,rand) result(x)
    real(defReal), intent(in)    :: P1
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: x
    real(defReal)                :: P1_loc
    real(defReal)                :: threshold
    integer(shortInt)            :: Low, Top, exec
    integer(shortInt), parameter :: UNIFORM = 1, LIN = 2, DELTA = 3
    character(100), parameter :: Here = 'sampleLegendre_P1 ( legendrePoly_func.f90)'

    ! Make local copy of P1 coeff. Take abs() to simplify code
    ! -ve P1 will be inverted at the end.
    P1_loc = abs(P1)

    ! Depending on whether P1 > 1 determine treshold and associated PDF for the mixing method
    ! For further details refer to Lux and Koblinger APPENDIX 3D
    ! If random number < threshold then Top is used.
    if ( P1_loc < ONE) then
      threshold = P1_loc
      Top = LIN
      Low = UNIFORM

    else if( P1_loc <= 3.0_defReal) then
      threshold = 0.5 * (P1_loc - ONE)
      Top = DELTA
      Low = LIN

    else
      call fatalError(Here,'P1 must have absolute value < 3.0')
      ! Avoid warnings
      threshold = ONE
      Top = 0
      Low = 0

    end if

    ! Use mixing method with the calculated Threshold
    if ( rand % get() < threshold ) then
      exec = Top
    else
      exec = Low
    end if

    ! Sample from UNIFORM ( PDF = 0.5); LIN ( PDF = 0.5 + 0.5 *mu) or DELTA ( PDF = DELTA(mu-1))
    select case(exec)
      case (UNIFORM)
        x = TWO * rand % get() - ONE

      case (LIN)
        ! Need to solve CDF(x) = 0.25 * x^2 + 0.5 * x + 0.25 = (0.5*x+0.5)^2)
        x = TWO * sqrt(rand % get()) - ONE

      case (DELTA)
        x = ONE

      case default
        call fatalError(Here,'This should never happen. WTF?')
        x = ZERO

    end select

    ! Invert result if P1 is -ve
    if ( P1 < ZERO ) x = -x

  end function sampleLegendre_P1

    
end module legendrePoly_func
