
add_clp_application(ImageInfo ImageInfo ${ITK_LIBRARIES})
add_clp_application(PrintImage PrintImage ${ITK_LIBRARIES} ${BLAS_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(4DImageMath 4DImageMath ${ITK_LIBRARIES} ${BLAS_LIBRARIES} ${EXPRTK_LIBRARY})
add_clp_application(VectorImageNormalization VectorImageNormalization ${ITK_LIBRARIES})
add_clp_application(ImageMultiplication ImageMultiplication ${ITK_LIBRARIES})

add_clp_application(4DToVectorImageConverter 4DToVectorImageConverter ${ITK_LIBRARIES})
add_clp_application(VectorTo4DImageConverter VectorTo4DImageConverter ${ITK_LIBRARIES})

add_clp_application(TextFileOperator TextFileOperator)