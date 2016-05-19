# manually trigger build of DD4Hep plugin
option (BUILD_DD4HEP_PLUGIN "build DD4HEP plugins" OFF)

if (BUILD_DD4HEP_PLUGIN)
  # we need DD4Hep
  find_package (DD4hep QUIET COMPONENTS DDCore DDSegmentation)

  if (DD4hep_FOUND)
    # get list of all source files
    file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

    # define library target
    add_library (ACTSDD4hepPlugin SHARED ${src_files})

    # setup include directories
    target_include_directories (ACTSDD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
    target_include_directories (ACTSDD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})

    # setup linker dependencies
    target_link_libraries (ACTSDD4hepPlugin PUBLIC ACTSCore)
    target_link_libraries (ACTSDD4hepPlugin PUBLIC ACTSTGeoPlugin)

    # set installation directories
    install (TARGETS ACTSDD4hepPlugin EXPORT ACTSTargets
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

    install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  # do not fail if ROOT is not available, just inform and continue
  else (DD4hep_FOUND)
    message (STATUS "${Red}DD4hep was not found -> skip build of DD4hepPlugins${ColorReset}")
  endif (DD4hep_FOUND)
else (BUILD_DD4HEP_PLUGIN)
  message (STATUS "${Blue}disable build of DD4hepPlugins${ColorReset}")
endif (BUILD_DD4HEP_PLUGIN)
