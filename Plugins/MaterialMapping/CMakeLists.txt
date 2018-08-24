file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (MaterialMappingPlugin SHARED ${src_files})

target_include_directories (MaterialMappingPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (MaterialMappingPlugin PUBLIC Acts::Core Acts::Legacy)

install (TARGETS MaterialMappingPlugin
         EXPORT ActsMaterialMappingTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT MaterialMappingPlugin TARGETS MaterialMappingPlugin)
