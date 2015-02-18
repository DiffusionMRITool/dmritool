
set(ExternalTAG_DIR ${DMRITOOL_DOXYGEN_OUTPUT_DIR})

find_package(PythonInterp REQUIRED)

add_custom_command( OUTPUT ${ExternalTAG_DIR}/InsightDoxygen.tag ${ExternalTAG_DIR}/vtkNightlyDoc.tag
  COMMAND ${CMAKE_COMMAND} -DITKTAG_DIR="${ExternalTAG_DIR}" -DVTKTAG_DIR="${ExternalTAG_DIR}" -P ${PROJECT_SOURCE_DIR}/CMake/DownloadDoxygenTAG.cmake
  COMMAND ${PYTHON_EXECUTABLE} "${PROJECT_SOURCE_DIR}/Utilities/GUnzip.py" "${ExternalTAG_DIR}/InsightDoxygen.tag.gz"
  COMMAND ${PYTHON_EXECUTABLE} "${PROJECT_SOURCE_DIR}/Utilities/GUnzip.py" "${ExternalTAG_DIR}/vtkNightlyDoc.tag.gz"
  COMMENT "Downloading and unpacking the Doxygen TAG"
  )
add_custom_target( ExternalDoxygenTAG ${CMAKE_COMMAND} -E echo "Finished obtaining external Doxygen TAG"
  DEPENDS 
  ${ExternalTAG_DIR}/InsightDoxygen.tag 
  ${ExternalTAG_DIR}/vtkNightlyDoc.tag
  )
