module hexagonalDeformationField_class

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
  type, public, extends(deformationField) :: hexagonalDeformationField
    
  contains
    ! Superclass procedures
    procedure :: forward
    procedure :: backward

  end type hexagonalDeformationField

contains

  !!
  !! Get value of the field at the co-ordinate point
  !!
  !!
  function forward(self, coords, X) result(r)
    class(hexagonalDeformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(:), intent(in), optional :: X
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta, hexRadius, hexRadiusOuter
    real(defReal)                     :: circumradius, deformedCircumradius, angle
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    
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

    if (present(X)) then
      delta = X(1) * self % delta(localID)
    else
      delta = self % delta(localID)
    end if

    circumradius = self % rS
    deformedCircumradius = self % rS + delta

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    angle = r(2)
    delta = sqrt(3.0) / TWO * (deformedCircumradius - circumradius) / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)

    ! For hexagonal perturbation, we check if within the hexagon
    hexRadius = sqrt(3.0_defReal) / TWO * self % rS / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)
    hexRadiusOuter = sqrt(3.0_defReal) / TWO * self % r0 / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)

    if (delta /= 0.0) then 
      if (r(1) < hexRadius) then 
        r(1) = r(1) + delta / hexRadius * r(1)
      elseif (r(1) > hexRadius .and. r(1) < hexRadiusOuter) then 
        r(1) = r(1) + delta / (hexRadius - hexRadiusOuter) * (r(1) - hexRadiusOuter)
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
    class(hexagonalDeformationField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal), dimension(:), intent(in), optional :: X
    real(defReal), dimension(3)       :: r
    real(defReal)                     :: delta, hexRadius, hexRadiusOuter, bound
    real(defReal)                     :: circumradius, deformedCircumradius, angle
    integer(shortInt)                 :: localID, temp, base
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: r_bar
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)

    ! Catch case if particle is outside the lattice
    ! Return unchanged position, as deformation is only applied within lattice
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

    if (present(X)) then
      delta = X(1) * self % delta(localID)
    else
      delta = self % delta(localID)
    end if

    circumradius = self % rS
    deformedCircumradius = self % rS + delta

    ! Transform to cylindrical coordinates, apply deformation, then transform back to Cartesian
    r = self % rTransform(r_bar)
    angle = r(2)
    delta = sqrt(3.0) / TWO * (deformedCircumradius - circumradius) / cos(modulo(angle, PI / 3.0) - PI / 6.0)
    ! For hexagonal perturbation, we check the radius in the hexagonal metric, i.e. the distance to the nearest hexagonal edge
    bound = sqrt(3.0_defReal) / TWO * deformedCircumradius / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)
    hexRadius = sqrt(3.0_defReal) / TWO * self % rS / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)
    hexRadiusOuter = sqrt(3.0_defReal) / TWO * self % r0 / cos(modulo(angle, PI / 3.0_defReal) - PI / 6.0_defReal)
    if (delta /= 0.0) then 
      if (r(1) < bound) then 
        r(1) = r(1) / (1 + delta / hexRadius)
      elseif (r(1) > bound .and. r(1) < hexRadiusOuter) then 
        r(1) = (r(1) + delta * hexRadiusOuter / (hexRadius - hexRadiusOuter)) / (1 + delta / (hexRadius - hexRadiusOuter))
      end if
    end if
        
    r = self % xTransform(r)

    ! Revert to original coordinates
    r = r + self % corner + (ijk - HALF) * self % pitch

  end function backward
    

end module hexagonalDeformationField_class
