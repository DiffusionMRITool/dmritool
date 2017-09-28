cmake_minimum_required(VERSION 2.8.8)
project(exprtk_builder C CXX)
include(ExternalProject)

ExternalProject_Add(exprtk
  GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git 
  CMAKE_ARGS ""
  CONFIGURE_COMMAND ""
  PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
  # Disable install step
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
  )

# Specify include dir
ExternalProject_Get_Property(exprtk source_dir)
message("source_dir=${source_dir}")
set(EXPRTK_INCLUDE_DIRS ${source_dir} PARENT_SCOPE)
set(EXPRTK_INCLUDE_DIRS ${source_dir})


message("EXPRTK_INCLUDE_DIRS=${EXPRTK_INCLUDE_DIRS}")
include_directories(${EXPRTK_INCLUDE_DIRS})

# add_subdirectory(test)

