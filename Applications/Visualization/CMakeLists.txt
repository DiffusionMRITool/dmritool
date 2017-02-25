
add_clp_application(OrientationsViewer OrientationsViewer ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_clp_application(MeshFromSHCoefficients MeshFromSHCoefficients ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(MeshFromSphericalFunctionTessellatedSamples MeshFromSphericalFunctionTessellatedSamples ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})