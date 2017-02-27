
add_clp_application(OrientationStatistics OrientationStatistics ${ITK_LIBRARIES} ${BLAS_LIBRARIES})
add_clp_application(SamplingSchemeQSpaceIncrementalEstimation SamplingSchemeQSpaceIncrementalEstimation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_clp_application(SamplingSchemeQSpaceIMOCEstimation SamplingSchemeQSpaceIMOCEstimation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_clp_application(SamplingSchemeQSpace1OptEstimation SamplingSchemeQSpace1OptEstimation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_clp_application(SamplingSchemeQSpaceRandomEstimation SamplingSchemeQSpaceRandomEstimation ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})

add_clp_application(SamplingSchemeDistance SamplingSchemeDistance ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
