
ADD_EXECUTABLE(test_spamsCppLasso cppLasso.cpp)
TARGET_LINK_LIBRARIES(test_spamsCppLasso ${ITK_LIBRARIES} ${BLAS_LIBRARIES})

add_test_application(test_spamsLasso spamsLassoTest ${ITK_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_test_application(test_spamsMatrix spamsMatrixTest ${ITK_LIBRARIES} ${BLAS_LIBRARIES})
