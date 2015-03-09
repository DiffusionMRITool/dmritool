
add_test_application(utlVNLBlasTest utlVNLBlasTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} )

add_gtest_application(utlDMRIGtest utlDMRIGTest ${BLAS_LIBRARIES})
add_gtest_application(utlCoreGtest utlCoreGTest)
add_gtest_application(utlCoreMKLGtest utlCoreMKLGTest ${MKL_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlGtest utlGTest ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_gtest_application(utlMathGtest utlMathGTest)
add_gtest_application(utlVNLBlasGTest utlVNLBlasGTest ${BLAS_LIBRARIES} ${ITK_LIBRARIES} )
add_gtest_application(utlVNLLapackGTest utlVNLLapackGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} )

add_gtest_application(utlNDArrayGTest utlNDArrayGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlVectorGTest utlVectorGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_gtest_application(utlMatrixGTest utlMatrixGTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})

# add_gtest_application(utlANNGTest utlANNGTest ANN ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
