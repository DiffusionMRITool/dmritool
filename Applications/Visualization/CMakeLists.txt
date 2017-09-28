
add_clp_application(OrientationsViewer OrientationsViewer ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_clp_application(MeshFromSHCoefficients MeshFromSHCoefficients ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(MeshFromDiscreteFiberODF MeshFromDiscreteFiberODF ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(MeshFromSphericalFunctionTessellatedSamples MeshFromSphericalFunctionTessellatedSamples ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(MeshFromPeaks MeshFromPeaks ${BLAS_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
# add_clp_application(VectorImageToVTK VectorImageToVTK ${ITK_LIBRARIES} ${BLAS_LIBRARIES})
add_clp_application(VectorImageToVTKXMLImageConverter VectorImageToVTKXMLImageConverter ${ITK_LIBRARIES} ${VTK_LIBRARIES})
add_clp_application(RGBImageToVTKXMLImageConverter RGBImageToVTKXMLImageConverter ${ITK_LIBRARIES} ${VTK_LIBRARIES})
add_clp_application(VTKXMLImageToVectorImageConverter VTKXMLImageToVectorImageConverter ${ITK_LIBRARIES} ${VTK_LIBRARIES})
