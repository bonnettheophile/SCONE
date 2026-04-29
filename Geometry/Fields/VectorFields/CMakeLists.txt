
add_sources( ./vectorField_inter.f90
             ./uniformVectorField_class.f90
             ./uniFissSitesField_class.f90
             ./weightWindowsField_class.f90
             ./deformationField_inter.f90
             ./radialDeformationField_class.f90
             ./hexagonalDeformationField_class.f90
             ./holeRadialDeformationField_class.f90
           )

add_unit_tests( ./Tests/uniformVectorField_test.f90
                ./Tests/uniFissSitesField_test.f90)

add_integration_tests( ./Tests/weightWindowsField_iTest.f90)
