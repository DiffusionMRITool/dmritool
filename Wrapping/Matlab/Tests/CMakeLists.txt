
add_executable(test_purecpp test_purecpp.cxx)
target_link_libraries(test_purecpp  ${ITK_LIBRARIES} ${GSL_LIBRARIES}) 

add_mex_file(test_puremex test_puremex )
add_mex_file(test_puremexitk test_puremexitk ${ITK_LIBRARIES})
add_mex_file(test_Spams test_Spams ${ITK_LIBRARIES} ${BLAS_LIBRARIES})
add_mex_file(test_Spams2 test_Spams2 ${ITK_LIBRARIES} ${BLAS_LIBRARIES})


if (DMRITOOL_USE_OPENMP)
  add_mex_file(test_puremexOpenmp test_puremexOpenmp )
endif()
