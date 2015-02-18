function(add_clp_test_application appName mainFile)
  set(${mainFile}_SOURCE ${mainFile}.cxx ${mainFile}.xml)
  GENERATECLP(${mainFile}_SOURCE ${mainFile}.xml "")

  add_executable(${appName} ${mainFile}.cxx)  
  target_link_libraries(${appName} ${ARGN})
endfunction(add_clp_test_application)


function(dmritool_install target)
  install(TARGETS ${target}
    RUNTIME DESTINATION dmritool/bin
    LIBRARY DESTINATION dmritool/lib
    ARCHIVE DESTINATION dmritool/lib/static
    BUNDLE DESTINATION /Applications/dmritool)
endfunction(dmritool_install)


function(add_clp_application appName mainFile)
  set(${mainFile}_SOURCE ${mainFile}.cxx)
  GENERATECLP(${mainFile}_SOURCE ${mainFile}.xml "")

  add_executable(${appName} ${mainFile}.cxx ${mainFile}.xml)
  target_link_libraries(${appName} ${ARGN})

  dmritool_install(${appName})
endfunction(add_clp_application)


function(add_application appName mainFile)
  add_executable(${appName} ${mainFile}.cxx)
  target_link_libraries(${appName} ${ARGN})
  dmritool_install(${appName})
endfunction(add_application)

function(add_test_application appName mainFile)
  add_executable(${appName} ${mainFile}.cxx)
  target_link_libraries(${appName} ${ARGN})
endfunction(add_test_application)

function(add_gtest_application appName mainFile)
  include_directories(${GTEST_INCLUDE_DIRS})
  add_executable(${appName} ${mainFile}.cxx)
  add_dependencies(${appName} googletest)

  if(NOT WIN32)
    target_link_libraries(${appName}
      ${GTEST_LIBS_DIR}/libgtest.a
      ${GTEST_LIBS_DIR}/libgtest_main.a
      pthread
      ${ARGN}
      )
  else()
    target_link_libraries(${appName}
      debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
      optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
      ${ARGN}
      )
    target_link_libraries(${appName}
      debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
      optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
      ${ARGN}
      )
  endif()

  if (BUILD_TESTING)
    add_test(${appName}_ctest ${appName})
  endif()
endfunction(add_gtest_application)


# copy folders
macro(configure_files srcDir destDir)
  if(${srcDir} IS_NEWER_THAN ${destDir})
    message(STATUS "Configuring directory ${destDir}")
    make_directory(${destDir})

    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/*)
    foreach(templateFile ${templateFiles})
      set(srcTemplatePath ${srcDir}/${templateFile})
      if(NOT IS_DIRECTORY ${srcTemplatePath})
        if(${srcTemplatePath} IS_NEWER_THAN ${destDir}/${templateFile})
          message(STATUS "Configuring file ${templateFile}")
          configure_file(${srcTemplatePath}  ${destDir}/${templateFile} COPYONLY)
        endif()
      else()
        configure_files(${srcTemplatePath} ${destDir}/${templateFile})
      endif()
    endforeach(templateFile)
  endif()
endmacro(configure_files)


# get all children folders
macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

