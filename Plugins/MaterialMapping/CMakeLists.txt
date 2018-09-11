file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsMaterialMappingPlugin SHARED ${src_files})

target_include_directories (ActsMaterialMappingPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ActsMaterialMappingPlugin PUBLIC ActsCore)

install (TARGETS ActsMaterialMappingPlugin
         EXPORT ActsMaterialMappingPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsMaterialMappingPlugin TARGETS ActsMaterialMappingPlugin)
