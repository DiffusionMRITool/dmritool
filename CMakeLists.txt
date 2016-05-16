cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

if(COMMAND cmake_policy)
  CMAKE_POLICY(VERSION 2.8.8)
endif(COMMAND cmake_policy)

project(DMRITOOL)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/CMake)
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(ProjectCMakeFunctions)

set(DMRITOOL_VERSION_MAJOR "0")
set(DMRITOOL_VERSION_MINOR "1")
set(DMRITOOL_VERSION_PATCH "1")

option(BUILD_STATIC_EXECUTABLES "Build static executables" OFF)
if(BUILD_STATIC_EXECUTABLES)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  set(BLA_STATIC TRUE)
else(BUILD_STATIC_EXECUTABLES)
  option(BUILD_SHARED_LIBS "Build shared libs" ON)
endif(BUILD_STATIC_EXECUTABLES)

set(DMRITOOL_SETTINGS_DIRECTORY "$ENV{HOME}/.dmritool")

if(NOT VERBOSITY_LEVEL)
  set(VERBOSITY_LEVEL 0 CACHE STRING
    "VERBOSITY Level for DEBUG."
    FORCE)
endif(NOT VERBOSITY_LEVEL)
add_definitions(-DUTL_VERBOSITY=${VERBOSITY_LEVEL})

# Support for C++11
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -stdlib=libc++")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")

option(BUILD_DOXYGEN_DOCUMENTATION "Build doxygen documentation" OFF)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(VTK_LIBRARIES_back ${VTK_LIBRARIES})

find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})
set(ITK_LIBRARIES_back ${ITK_LIBRARIES})

find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})

set(VTK_LIBRARIES ${VTK_LIBRARIES_back})
set(ITK_LIBRARIES ${ITK_LIBRARIES_back})

find_path(MKL_INCLUDE_DIR mkl_cblas.h 
  /opt/intel/mkl/include
  $ENV{HOME}/intel/mkl/include
  /usr/include
  /usr/local/include
  )
mark_as_advanced(MKL_INCLUDE_DIR)
if (IS_DIRECTORY ${MKL_INCLUDE_DIR})
  option(DMRITOOL_USE_MKL "Use MKL for numerical algebra" ON)
else()
  option(DMRITOOL_USE_MKL "Use MKL for numerical algebra" OFF)
endif()

if (DMRITOOL_USE_MKL)
  find_package(MKL REQUIRED)
  message("MKL_LIBRARIES=${MKL_LIBRARIES}")
  add_definitions(-DUTL_USE_MKL)
  set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
  set(BLAS_LIBRARIES ${MKL_LIBRARIES})
else()
  ## BLAS - required by utl::NDArray, spams
  find_package(OpenBLAS REQUIRED)

  ## LAPACK - required by utl::NDArray, spams
  find_package(LAPACK REQUIRED)
  set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
  set(BLAS_LIBRARIES ${OpenBLAS_LIBRARY} ${LAPACK_LIBRARIES})
endif()
set(LAPACK_FOUND TRUE)
set(BLAS_FOUND TRUE)
if (BUILD_STATIC_EXECUTABLES)
  set(BLAS_LIBRARIES "${BLAS_LIBRARIES};-lgfortran")
endif()
message("BLAS_LIBRARIES=${BLAS_LIBRARIES}")
message("LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")

# message("ITK_LIBRARIES=${ITK_LIBRARIES}")
# message("VTK_LIBRARIES=${VTK_LIBRARIES}")

option(DMRITOOL_USE_FASTLAPACK "Use fast version of LAPACK routines, it needs to build openblas manually, or use MKL." OFF)
if (DMRITOOL_USE_FASTLAPACK)
  set(DMRITOOL_USE_FASTLAPACK ON)
  add_definitions(-DUTL_USE_FASTLAPACK)
endif()


## for spams
add_definitions(-DUSE_BLAS_LIB=1)

## Required Packages
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
message("GSL_LIBRARIES=${GSL_LIBRARIES}")

## OpenMP
option(DMRITOOL_USE_OPENMP "Use double precision OpenMP if found" OFF)
if(DMRITOOL_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  add_definitions(-DUTL_USE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


# Copy Data
configure_files(${DMRITOOL_SOURCE_DIR}/Data ${DMRITOOL_SETTINGS_DIRECTORY}/Data)

# Configuration files with build settings
configure_file(${DMRITOOL_SOURCE_DIR}/CMake/DMRITOOLConfigure.h.in ${DMRITOOL_BINARY_DIR}/DMRITOOLConfigure.h)

## BUILD_TESTING
option(BUILD_TESTING "Build tests" ON)
if (BUILD_TESTING)
  set(buildalltests "alltests")
  add_custom_target(${buildalltests}) 
  enable_testing()
  include(CTest)
endif()

# Google Test
if (BUILD_TESTING)
  add_subdirectory(ThirdParty/gtest)
  message("GTEST_INCLUDE_DIRS=${GTEST_INCLUDE_DIRS}")
  message("GTEST_LIBS_DIR=${GTEST_LIBS_DIR}")
  include_directories(${GTEST_INCLUDE_DIRS})
endif ()


## Modules
subdirlist(DMRITOOL_MODULE_DIRECTORIES ${DMRITOOL_SOURCE_DIR}/Modules)


set(THIRDPARTY_DIRECTORIES
  Spams
  )


set(DMRITOOL_INCLUDE_DIRECTORIES ${DMRITOOL_BINARY_DIR} ${DMRITOOL_SOURCE_DIR})

foreach(dir ${DMRITOOL_MODULE_DIRECTORIES})
  set(DMRITOOL_INCLUDE_DIRECTORIES ${DMRITOOL_INCLUDE_DIRECTORIES} ${DMRITOOL_SOURCE_DIR}/Modules/${dir}/include)
endforeach(dir)


foreach(dir ${THIRDPARTY_DIRECTORIES})
  set(DMRITOOL_INCLUDE_DIRECTORIES ${DMRITOOL_INCLUDE_DIRECTORIES} ${DMRITOOL_SOURCE_DIR}/ThirdParty/${dir}/include)
endforeach(dir)

include_directories(${DMRITOOL_INCLUDE_DIRECTORIES})


foreach(dir ${DMRITOOL_MODULE_DIRECTORIES})
  set(dir_abs ${PROJECT_SOURCE_DIR}/Modules/${dir})
  if (EXISTS "${dir_abs}}/src" AND EXISTS "${dir_abs}/src/CMakeLists.txt")
    add_subdirectory(${dir_abs}/src)
  endif()
  if (BUILD_TESTING AND EXISTS "${dir_abs}/test" AND EXISTS "${dir_abs}/test/CMakeLists.txt")
    add_subdirectory(${dir_abs}/test)
  endif()
endforeach(dir)

foreach(dir ${THIRDPARTY_DIRECTORIES})
  set(dir_abs ${PROJECT_SOURCE_DIR}/ThirdParty/${dir})
  if (EXISTS "${dir_abs}}/src" AND EXISTS "${dir_abs}/src/CMakeLists.txt")
    add_subdirectory(${dir_abs}/src)
  endif()
  if(BUILD_TESTING AND EXISTS "${dir_abs}/test" AND EXISTS "${dir_abs}/test/CMakeLists.txt")
    add_subdirectory(${dir_abs}/test)
  endif ()
endforeach(dir)


##
# Ensure inclusion of files in Xcode/Qt project file.
#
if(APPLE)
  file(GLOB_RECURSE HEADER_FILES ${DMRITOOL_SOURCE_DIR}/Modules/*.h)
  file(GLOB_RECURSE IMPLEMENTATION_FILES ${DMRITOOL_SOURCE_DIR}/Modules/*.hxx)

  add_library(MODULES MODULE ${HEADER_FILES} ${IMPLEMENTATION_FILES})
  target_link_libraries(MODULES ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES})
endif(APPLE)


##
# Applications
#
add_subdirectory(Applications)
add_subdirectory(Utilities)

if (BUILD_QT_APPLICATIONS)
  find_package(Qt4 REQUIRED)
else()
  find_package(Qt4)
endif()
if (NOT QT4_FOUND AND NOT QT_FOUND)
  option(BUILD_QT_APPLICATIONS "Build QT Applications" OFF)
else()
  option(BUILD_QT_APPLICATIONS "Build QT Applications" ON)
  # Widgets finds its own dependencies.
  # find_package(Qt5Widgets REQUIRED)
  set(QT_USE_QTSVG TRUE)
  set(QT_USE_QTXML TRUE)

  include(${QT_USE_FILE})

  if (BUILD_QT_APPLICATIONS)
    add_subdirectory(QTApplications)
  endif()
endif()


##
# Wrapping 
#
add_subdirectory(Wrapping)

## Testing
IF (BUILD_TESTING)
  add_subdirectory(Testing)
endif()

##
# Output Directories for Executable, Library, and Archive Targets
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${DMRITOOL_BINARY_DIR}/bin CACHE INTERNAL 
  "Executable output directory")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${DMRITOOL_BINARY_DIR}/lib CACHE INTERNAL
  "Library output directory")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${DMRITOOL_BINARY_DIR}/lib CACHE INTERNAL
  "Archive output directory")

set(DMRITOOL_LIBRARY_DIRECTORIES ${DMRITOOL_BINARY_DIR}/lib)
set(DMRITOOL_USE_FILE ${DMRITOOL_BINARY_DIR}/UseDMRITOOL.cmake)
configure_file(CMake/DMRITOOLConfig.cmake.in ${DMRITOOL_BINARY_DIR}/DMRITOOLConfig.cmake)
configure_file(CMake/UseDMRITOOL.cmake.in  ${DMRITOOL_BINARY_DIR}/UseDMRITOOL.cmake COPYONLY IMMEDIATE)

##
# Doxygen Documentation
#
if(BUILD_DOXYGEN_DOCUMENTATION)
    
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DMRITOOL_DOXYGEN_OUTPUT_DIR ${DMRITOOL_BINARY_DIR}/Documentation/Doxygen)
    configure_file(${DMRITOOL_SOURCE_DIR}/Documentation/Doxygen/doxygen.in ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/doxygen.config)

    option(LINK_EXTERNAL_DOXYGEN "Build doxygen document with external tags" OFF)
    if(LINK_EXTERNAL_DOXYGEN)
      include(ExternalDoxygenTAG)
      add_custom_target(doxygen 
        COMMAND sed -e 's/PROJECT_NUMBER.*=.*\\>/PROJECT_NUMBER = '$$\( cd ${DMRITOOL_SOURCE_DIR} && git describe --abbrev=8 --always\)'/' ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/doxygen.config | doxygen -
        # COMMAND ${DOXYGEN} ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/doxygen.config
        DEPENDS ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/InsightDoxygen.tag ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/vtkNightlyDoc.tag
        )
    else()
      add_custom_target(doxygen 
        COMMAND sed -e 's/PROJECT_NUMBER.*=.*\\>/PROJECT_NUMBER = '$$\( cd ${DMRITOOL_SOURCE_DIR} && git describe --abbrev=8 --always\)'/' ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/doxygen.config | doxygen -
        # COMMAND ${DOXYGEN} ${DMRITOOL_DOXYGEN_OUTPUT_DIR}/doxygen.config
        )
    endif()
  endif(DOXYGEN_FOUND)

endif(BUILD_DOXYGEN_DOCUMENTATION)


##
# Packaging
#
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_README  "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${DMRITOOL_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${DMRITOOL_VERSION_MINOR}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DMRITOOL")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "dmritool-${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

include(CPack)


