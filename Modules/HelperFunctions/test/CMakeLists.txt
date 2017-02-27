
add_test_application(utlVNLBlasTest utlVNLBlasTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} )
add_test_application(utlGradientTest utlGradientTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} )

add_gtest_application(utlDMRIGTest utlDMRIGTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlCoreGTest utlCoreGTest)
add_gtest_application(utlCoreMKLGTest utlCoreMKLGTest ${MKL_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlGTest utlGTest ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_gtest_application(utlMathGTest utlMathGTest)
add_gtest_application(utlVNLBlasGTest utlVNLBlasGTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} )
add_gtest_application(utlVNLLapackGTest utlVNLLapackGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} )

add_gtest_application(utlNDArrayGTest utlNDArrayGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlVectorGTest utlVectorGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlMatrixGTest utlMatrixGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utl4thOrderTensorGTest utl4thOrderTensorGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})

# add_gtest_application(utlANNGTest utlANNGTest ANN ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
