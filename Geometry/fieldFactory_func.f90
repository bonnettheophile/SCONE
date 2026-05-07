!!
!! Module to build new fields and add them to the geometry registry
!!
module fieldFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Fields
  use field_inter,                     only : field
  use uniformScalarField_class,        only : uniformScalarField
  use uniformVectorField_class,        only : uniformVectorField
  use uniFissSitesField_class,         only : uniFissSitesField
  use weightWindowsField_class,        only : weightWindowsField
  use cartesianField_class,            only : cartesianField
  use deformationField_inter,          only : deformationField
  use radialDeformationField_class,    only : radialDeformationField
  use hexagonalDeformationField_class, only : hexagonalDeformationField
  use holeRadialDeformationField_class, only : holeRadialDeformationField
  use axialDeformationField_class,      only : axialDeformationField
  use gapDeformationField_class,       only : gapDeformationField

  ! Geometry
  use geometryReg_mod,          only : gr_addField => addField

  implicit none
  private


  !! Parameters
  character(nameLen), dimension(*), parameter :: AVAILABLE_FIELDS = ['uniformScalarField        ',&
                                                                     'uniformVectorField        ',&
                                                                     'uniFissSitesField         ',&
                                                                     'weightWindowsField        ',&
                                                                     'cartesianField            ',&
                                                                     'radialDeformationField    ',&
                                                                     'hexagonalDeformationField ',&
                                                                     'holeRadialDeformationField',&
                                                                     'axialDeformationField     ',&
                                                                     'gapDeformationField       ']

   ! Public interface
   public :: new_field

contains

  !!
  !! Initialises field from dictionary. The field gets then allocated
  !! inside the geometry registry
  !!
  !! Args:
  !!   dict [in]   -> Dictionary with definition
  !!   name [in]   -> Name of the field for the geometry registry
  !!
  !! Errors:
  !!   fatalError if type of field is unknown
  !!
  subroutine new_field(dict, name)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    class(field), allocatable      :: kentta
    character(nameLen)             :: type
    character(100), parameter :: Here = 'new_field (fieldFactory_func.f90) '

    ! Get type
    call dict % get(type, 'type')

    ! Build Field
    select case (type)
      case ('uniformScalarField')
        allocate(uniformScalarField :: kentta)

      case ('uniformVectorField')
        allocate(uniformVectorField :: kentta)

      case ('uniFissSitesField')
        allocate(uniFissSitesField :: kentta)

      case ('weightWindowsField')
        allocate(weightWindowsField :: kentta)

      case ('cartesianField')
        allocate(cartesianField :: kentta)
      
      case ('radialDeformationField')
        allocate(radialDeformationField :: kentta)
      
      case ('hexagonalDeformationField')
        allocate(hexagonalDeformationField :: kentta)

      case ('holeRadialDeformationField')
        allocate(holeRadialDeformationField :: kentta)

      case ('axialDeformationField')
        allocate(axialDeformationField :: kentta)

      case ('gapDeformationField')
        allocate(gapDeformationField :: kentta)

      case default
        print '(A)', "AVAILABLE FIELDS:"
        print '(A)', AVAILABLE_FIELDS
        call fatalError(Here, trim(type)//' is not valid field. See list above.')

    end select

    ! Initialise Field
    call kentta % init(dict)

    ! Call geometry registry to add field
    call gr_addField(kentta, name)

  end subroutine new_field
  
end module fieldFactory_func
