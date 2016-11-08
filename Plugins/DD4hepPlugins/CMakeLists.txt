# manually trigger build of DD4Hep plugin
option (BUILD_DD4HEP_PLUGIN "build DD4HEP plugins" OFF)

if (BUILD_DD4HEP_PLUGIN)
  # we need the TGeoPlugin
  set (BUILD_TGEO_PLUGIN ON PARENT_SCOPE)

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
    install (TARGETS ACTSDD4hepPlugin EXPORT ACTSDD4hepPluginTargets
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

    install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    # add to list of supported components
    list (APPEND _supported_components DD4hepPlugin)
    set (_supported_components "${_supported_components}" PARENT_SCOPE)
    
    # add dummy project capturing all targets for the CDash DD4hepPlugin subproject
    add_custom_target(DD4hepPluginProject)
    add_dependencies(DD4hepPluginProject ACTSDD4hepPlugin)
    
  # do not fail if ROOT is not available, just inform and continue
  else (DD4hep_FOUND)
    message (STATUS "${Red}DD4hep was not found -> skip build of DD4hepPlugins${ColorReset}")
  endif (DD4hep_FOUND)
else (BUILD_DD4HEP_PLUGIN)
  message (STATUS "${Blue}disable build of DD4hepPlugins${ColorReset}")
endif (BUILD_DD4HEP_PLUGIN)
