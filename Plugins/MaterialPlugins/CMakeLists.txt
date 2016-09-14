# manually trigger build of Material plugin
option (BUILD_MATERIAL_PLUGIN "build Material plugins" OFF)

if (BUILD_MATERIAL_PLUGIN)

  # get list of all source files
  file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")
  # define library target
  add_library (ACTSMaterialPlugin SHARED ${src_files})

  # setup include directories
  target_include_directories (ACTSMaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

  # setup linker dependencies
  target_link_libraries (ACTSMaterialPlugin PUBLIC ACTSCore)

  target_compile_options (ACTSMaterialPlugin PUBLIC -std=c++14)
  # set installation directories
  install (TARGETS ACTSMaterialPlugin EXPORT ACTSMaterialPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

  install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  # add to list of supported components
  list (APPEND _supported_components MaterialPlugin)
  set (_supported_components "${_supported_components}" PARENT_SCOPE)
   
else (BUILD_MATERIAL_PLUGIN)
    message (STATUS "${Blue}disable build of MaterialPlugins${ColorReset}")
endif (BUILD_MATERIAL_PLUGIN)
