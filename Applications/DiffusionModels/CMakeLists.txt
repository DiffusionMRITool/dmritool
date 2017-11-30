
# Spherical Tessellation
add_clp_application(SphereTessellator SphereTessellator ${ITK_LIBRARIES})

#ADD_CUSTOM_COMMAND (TARGET SphereTessellator 
#  POST_BUILD
#  COMMAND mkdir -p ${DMRITOOL_BINARY_DIR}/Directions
#  COMMAND SphereTessellator --tessorder 1 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t1.txt
#  COMMAND SphereTessellator --tessorder 2 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t2.txt
#  COMMAND SphereTessellator --tessorder 3 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t3.txt
#  COMMAND SphereTessellator --tessorder 4 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t4.txt
#  COMMAND SphereTessellator --tessorder 5 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t5.txt
#  COMMAND SphereTessellator --tessorder 6 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t6.txt
#  COMMAND SphereTessellator --tessorder 7 --hemisphere ${DMRITOOL_BINARY_DIR}/Directions/directions_t7.txt
#  COMMENT "Generating gradient files in ${DMRITOOL_BINARY_DIR}/Directions/"
#  )

add_clp_application(SHCoefficientsRotate SHCoefficientsRotate ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SHCoefficientsPower SHCoefficientsPower ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SHCoefficientsToSphericalFunctionSamples SHCoefficientsToSphericalFunctionSamples ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SHCoefficientsToGFA SHCoefficientsToGFA ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES} )
add_clp_application(ODFNormalization ODFNormalization ${ITK_LIBRARIES} ${GSL_LIBRARIES})

add_clp_application(SphericalPolarFourierImaging SphericalPolarFourierImaging ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES} )
add_clp_application(MeanDiffusivityEstimator MeanDiffusivityEstimator ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SPFToProfile SPFToProfile ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SPFToODF SPFToODF ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(SPFToScalarMap SPFToScalarMap ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})

add_clp_application(ComputeSHCoefficientsOfDWIFromSymmetricTensor ComputeSHCoefficientsOfDWIFromSymmetricTensor ${ITK_LIBRARIES} ${GSL_LIBRARIES})

add_clp_application(TensorFormatConverter TensorFormatConverter ${BLAS_LIBRARIES} ${ITK_LIBRARIES})
