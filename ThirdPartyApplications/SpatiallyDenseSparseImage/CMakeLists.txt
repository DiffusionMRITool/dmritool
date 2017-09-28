PROJECT(SparseImage)

add_clp_application(SparseToVectorImageConverter SparseToVectorImageConverter ${ITK_LIBRARIES})
add_clp_application(VectorToSparseImageConverter VectorToSparseImageConverter ${ITK_LIBRARIES})


