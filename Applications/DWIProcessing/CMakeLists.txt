
add_clp_application(DWISimulator DWISimulator ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(DWISingleVoxelSimulator DWISingleVoxelSimulator ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES} ${GSL_LIBRARIES})
add_clp_application(DWINoiseGenerator DWINoiseGenerator ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})
add_clp_application(DWIMultipleShellSplit DWIMultipleShellSplit ${ITK_LIBRARIES})
add_clp_application(CorrectBValues CorrectBValues )
add_clp_application(DWIPreprocess DWIPreprocess ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${ITK_LIBRARIES})

