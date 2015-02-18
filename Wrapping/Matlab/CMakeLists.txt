project(DMRITOOL_MEX)

cmake_minimum_required(VERSION 2.8)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL 
  "Executable output directory")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL
  "Library output directory")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
include(ProjectCMakeFunctions)
include(MatlabMakeMacros)

option(BUILD_SHARED_LIBS "Build shared libs" ON)


find_package(Matlab REQUIRED)
include_directories(${MATLAB_INCLUDE_DIR})
message("MATLAB_LIBRARIES=${MATLAB_LIBRARIES}")

find_package(DMRITOOL REQUIRED)
include(${DMRITOOL_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if(DMRITOOL_USE_MKL)
  find_package(MKL REQUIRED)
  add_definitions(-DUTL_USE_MKL)
  set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
  set(BLAS_LIBRARIES ${MKL_LIBRARIES})
else()
  ## BLAS - required by spams
  find_package(OpenBLAS REQUIRED)
  set(BLAS_LIBRARIES ${OpenBLAS_LIBRARY})

  ## LAPACK - required by spams
  find_package(LAPACK REQUIRED)
  set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${OpenBLAS_LIBRARY};-llapacke")
  set(BLAS_LIBRARIES ${OpenBLAS_LIBRARY} ${LAPACK_LIBRARIES})
endif()
set(LAPACK_FOUND TRUE)
set(BLAS_FOUND TRUE)
add_definitions(-DUSE_BLAS_LIB=1)

if (DMRITOOL_USE_FASTLAPACK)
  add_definitions(-DUTL_USE_FASTLAPACK)
endif()

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

if(NOT VERBOSITY_LEVEL)
  set(VERBOSITY_LEVEL 0 CACHE STRING
    "VERBOSITY Level for DEBUG."
    FORCE)
endif(NOT VERBOSITY_LEVEL)
add_definitions(-DUTL_VERBOSITY=${VERBOSITY_LEVEL})

if(DMRITOOL_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  add_definitions(-DUTL_USE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



include_directories(${DMRITOOL_MEX_SOURCE_DIR}/MEXHelperFunctions)

add_subdirectory(Common)
add_subdirectory(DWIProcessing)
add_subdirectory(DiffusionModels)
add_subdirectory(SamplingScheme)

add_subdirectory(Tests)


