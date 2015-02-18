
FUNCTION(add_matlab_cpp_executable appName mainFile)

  # setup doc for matlab functions
  set(MATLAB_FUNCTION_NAME ${appName})
  configure_file(${appName}.m ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${appName}.m)

  add_library(${appName} ${mainFile} )
  target_link_libraries(${appName} ${MATLAB_LIBRARIES} ${ARGN})

  if(MSYS)
    target_link_libraries(${appName} stdc++)
  endif()

  add_custom_command(TARGET ${appName} 
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/bin
    COMMAND ${PROJECT_SOURCE_DIR}/utilsMEX/mex_link.sh  ${appName} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${PROJECT_SOURCE_DIR}/utilsMEX/mex_stub.cpp
    )

ENDFUNCTION(add_matlab_cpp_executable)


FUNCTION(add_matlab_c_executable appName mainFile)

  # setup doc for matlab functions
  set(MATLAB_FUNCTION_NAME ${appName})
  configure_file(${appName}.m ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${appName}.m)

  add_library(${appName} ${mainFile} )
  target_link_libraries(${appName} ${MATLAB_LIBRARIES} ${ARGN})

  if(MSYS)
    target_link_libraries(${appName} stdc++)
  endif()

  add_custom_command(TARGET ${appName} 
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/bin
    COMMAND ${PROJECT_SOURCE_DIR}/utilsMEX/mex_link.sh  ${appName} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${PROJECT_SOURCE_DIR}/utilsMEX/mex_stub.c
    )

ENDFUNCTION(add_matlab_c_executable)

