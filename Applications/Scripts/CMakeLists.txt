
file(GLOB pyfiles *.py)
file(GLOB shfiles *.sh)

set(copy_files ${pyfiles} ${shfiles})
# message("copy_files=${copy_files}")

file(MAKE_DIRECTORY ${DMRITOOL_BINARY_DIR}/bin)
file(COPY ${copy_files} DESTINATION ${DMRITOOL_BINARY_DIR}/bin)


