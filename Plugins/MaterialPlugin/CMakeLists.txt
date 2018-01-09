file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ACTSMaterialPlugin SHARED ${src_files})

target_include_directories (ACTSMaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ACTSMaterialPlugin PUBLIC ACTSCore)

install (TARGETS ACTSMaterialPlugin
         EXPORT ACTSMaterialPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT MaterialPlugin TARGETS ACTSMaterialPlugin)
