
add_test_application(itkUnaryFunctorLookUpTableTest itkUnaryFunctorLookUpTableTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES})
add_test_application(itkFunctorHashTableTest itkFunctorHashTableTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES} )

add_test_application(itkSphericalPolarFourierGeneratorTest itkSphericalPolarFourierGeneratorTest ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_test_application(itkSphericalHarmonicsGeneratorTest itkSphericalHarmonicsGeneratorTest ${ITK_LIBRARIES} ${GSL_LIBRARIES})

add_gtest_application(itkSphericalHarmonicsGeneratorGTest itkSphericalHarmonicsGeneratorGTest ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_gtest_application(itkUnaryFunctorLookUpTableGTest itkUnaryFunctorLookUpTableGTest ${ITK_LIBRARIES})
