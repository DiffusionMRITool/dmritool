
add_mex_file(mexGetSHBasisMatrix mexGetSHBasisMatrix.cxx ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_mex_file(mexGetSPFBasisMatrix mexGetSPFBasisMatrix.cxx ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})

add_mex_file(mexSphericalPolarFourierImaging mexSphericalPolarFourierImaging.cxx ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES} )
add_mex_file(mexGetSamplesFromSPFCoefficients mexGetSamplesFromSPFCoefficients.cxx  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
