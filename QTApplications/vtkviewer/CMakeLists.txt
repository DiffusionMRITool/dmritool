# VTK Viewer
# Written 2012 Hal Canary <http://cs.unc.edu/~hal>
# Copyright 2012 University of North Carolina at Chapel Hill.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you
# may not use this file except in compliance with the License.  You
# may obtain a copy of the License at
#
#   LICENSE.md in this repository or
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.  See the License for the specific language governing
# permissions and limitations under the License.
cmake_minimum_required(VERSION 2.6)
project(vtkviewer)
# set(CMAKE_BUILD_TYPE debug)
if (NOT VTK_FOUND)
  set(VTK_DIR "VTK_DIR-NOTFOUND" CACHE PATH "location of VTK libraries" )
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()
set(QT_USE_OPENGL TRUE)
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
qt4_wrap_cpp(MOC_GENERATED VTKViewer.h)
add_executable(vtkviewer main.cxx VTKViewer.cxx ${MOC_GENERATED})
target_link_libraries(vtkviewer ${VTK_LIBRARIES} ${QT_QTOPENGL_LIBRARY})
