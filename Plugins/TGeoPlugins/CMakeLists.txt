# manually trigger build of TGeo plugin
option (BUILD_TGEO_PLUGIN "build TGeo plugins" OFF)

if (BUILD_TGEO_PLUGIN)
  # we need root
  find_package (ROOT QUIET COMPONENTS Geom)

  if (ROOT_FOUND)
    # get list of all source files
    file (GLOB src_files "src/*.cpp")

    # define library target
    add_library (ACTSTGeoPlugin SHARED ${src_files})

    # setup include directories
    target_include_directories (ACTSTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
    target_include_directories (ACTSTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

    # setup linker dependencies
    target_link_libraries (ACTSTGeoPlugin PUBLIC ACTSCore)
    target_link_libraries (ACTSTGeoPlugin PUBLIC Geom)

    # set installation directories
    install (TARGETS ACTSTGeoPlugin EXPORT ACTSTargets
             LIBRARY DESTINATION lib)

    install (DIRECTORY include/ACTS DESTINATION include)
  # do not fail if ROOT is not available, just inform and continue
  else (ROOT_FOUND)
    message (STATUS "ROOT was not found -> skip build of TGeoPlugins")
  endif (ROOT_FOUND)
endif (BUILD_TGEO_PLUGIN)
