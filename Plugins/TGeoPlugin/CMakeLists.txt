if (BUILD_TGEO_PLUGIN)
  # we need root
  find_package (ROOT ${ACTS_ROOT_VERSION} ${ACTS_FIND_QUIET_TEXT} COMPONENTS Geom)

  if (ROOT_FOUND)
    # get list of all source files
    file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

    # define library target
    add_library (ACTSTGeoPlugin SHARED ${src_files})

    # setup include directories
    target_include_directories (ACTSTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
    target_include_directories (ACTSTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

    # setup linker dependencies
    # temporary support for ROOT version 6.08 without CMake Namespace
    target_link_libraries (ACTSTGeoPlugin PUBLIC ACTSCore)
    if (${ROOT_VERSION_MINOR} MATCHES "8" )
      target_link_libraries (ACTSTGeoPlugin PUBLIC Geom)
    else (${ROOT_VERSION_MINOR} MATCHES "8" )
      target_link_libraries (ACTSTGeoPlugin PUBLIC ROOT::Geom)
    endif (${ROOT_VERSION_MINOR} MATCHES "8" )

    # set installation directories
    install (TARGETS ACTSTGeoPlugin EXPORT ACTSTGeoPluginTargets
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

    install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    # add to list of supported components
    list (APPEND _supported_components TGeoPlugin)
    set (_supported_components "${_supported_components}" PARENT_SCOPE)
    
    # add to CDash TGeoPlugin subproject
    acts_add_targets_to_cdash_project(PROJECT TGeoPlugin TARGETS ACTSTGeoPlugin)
  # do not fail if ROOT is not available, just inform and continue
  else (ROOT_FOUND)
    message (STATUS "${Red}ROOT was not found -> skip build of TGeoPlugins${ColorReset}")
  endif (ROOT_FOUND)
else (BUILD_TGEO_PLUGIN)
  message (STATUS "${Blue}disable build of TGeoPlugins${ColorReset}")
endif (BUILD_TGEO_PLUGIN)
