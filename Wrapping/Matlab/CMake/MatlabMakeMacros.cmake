# This file was found in the source code of Dru F., Fillard P.,
# Vercauteren T., "An ITK Implementation of the Symmetric Log-Domain
# Diffeomorphic Demons Algorithm", Insight Journal, 2009 Jan-Jun
# http://hdl.handle.net/10380/3060
#
# The file has been edited for project Gerardus
# Version: 0.2.0
# $Rev$
# $Date$

MACRO(ADD_MEX_FILE Target mainFile)
  INCLUDE_DIRECTORIES("${MATLAB_INCLUDE_DIR}")

  ADD_LIBRARY(${Target} SHARED ${mainFile})
  set_target_properties(${Target} PROPERTIES COMPILE_FLAGS "-DMATLAB_MEX_FILE")
  
  IF(WIN32)
    TARGET_LINK_LIBRARIES(${Target} 
      ${MATLAB_LIBRARIES} 
      )
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES(${Target} 
      ${MATLAB_LIBRARIES} 
      m
      )
  ENDIF(WIN32)
    
  TARGET_LINK_LIBRARIES(${Target} ${ARGN}) 
  
  SET_TARGET_PROPERTIES(${Target} PROPERTIES PREFIX "")
  
  # Determine mex suffix
  IF(UNIX)
    # if this is OSX (which is UNIX) then the library suffixes depend on the architecture
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac intel 32-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmaci")
      ELSEIF(CMAKE_OSX_ARCHITECTURES MATCHES x86_64)
        # mac intel 64-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmaci64")
      ELSE(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # Mac Power PC
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmac")
      ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexglx")
      ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexa64")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
        MESSAGE(FATAL_ERROR 
          "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ELSEIF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexw32")
    ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
      SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexw64")
    ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
      MESSAGE(FATAL_ERROR 
        "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
    ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
  ENDIF(UNIX)
  
  IF(MSVC)
    SET(MATLAB_FLAGS "-DMATLAB_MEX_FILE")
    SD_APPEND_TARGET_PROPERTIES(${Target} COMPILE_FLAGS "${MATLAB_FLAGS}")
    SET_TARGET_PROPERTIES(${Target} PROPERTIES LINK_FLAGS "/export:mexFunction")
  ELSE(MSVC)
    IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET(MATLAB_FLAGS "-fPIC" "-D_GNU_SOURCE" "-pthread"
        "-D_FILE_OFFSET_BITS=64" "-DMX_COMPAT_32")
    ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET(MATLAB_FLAGS "-fPIC" "-D_GNU_SOURCE" "-pthread"
        "-D_FILE_OFFSET_BITS=64")
    ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    SD_APPEND_TARGET_PROPERTIES(${Target} COMPILE_FLAGS "${MATLAB_FLAGS}")
    
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac intel 32-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_ROOT}/bin/maci -Wl,-flat_namespace -undefined suppress")
      ELSEIF(CMAKE_OSX_ARCHITECTURES MATCHES x86_64)
        # mac intel 64-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_ROOT}/bin/maci64 -Wl,-flat_namespace -undefined suppress")
      ELSE(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac powerpc?
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_SYS} -Wl,-flat_namespace -undefined suppress")
      ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-Wl,-E -Wl,--no-undefined")
      ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-Wl,-E -Wl,--no-undefined")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
        MESSAGE(FATAL_ERROR 
          "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ENDIF(MSVC)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${Target}.m")
    set(MATLAB_FUNCTION_NAME ${Target})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${Target}.m" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${Target}.m)
  endif()

ENDMACRO(ADD_MEX_FILE)


MACRO(SD_APPEND_TARGET_PROPERTIES TARGET_TO_CHANGE PROP_TO_CHANGE)
  FOREACH(_newProp ${ARGN})
    GET_TARGET_PROPERTY(_oldProps ${TARGET_TO_CHANGE} ${PROP_TO_CHANGE})
    IF(_oldProps)
      IF(NOT "${_oldProps}" MATCHES "^.*${_newProp}.*$")
        SET_TARGET_PROPERTIES(${TARGET_TO_CHANGE} PROPERTIES ${PROP_TO_CHANGE} "${_newProp} ${_oldProps}")
      ENDIF(NOT "${_oldProps}" MATCHES "^.*${_newProp}.*$")
    ELSE(_oldProps)
      SET_TARGET_PROPERTIES(${TARGET_TO_CHANGE} PROPERTIES ${PROP_TO_CHANGE} ${_newProp})
    ENDIF(_oldProps)
  ENDFOREACH(_newProp ${ARGN})
ENDMACRO(SD_APPEND_TARGET_PROPERTIES TARGET_TO_CHANGE PROP_TO_CHANGE)



