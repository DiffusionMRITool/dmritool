
set(copy_files 
  VTKPolyData.py
  )

file(MAKE_DIRECTORY ${DMRITOOL_BINARY_DIR}/bin)
file(COPY ${copy_files} DESTINATION ${DMRITOOL_BINARY_DIR}/bin)
