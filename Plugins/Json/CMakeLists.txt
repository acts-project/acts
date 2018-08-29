file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsJsonPlugin SHARED ${src_files})

target_include_directories (ActsJsonPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ActsJsonPlugin PUBLIC ActsCore)

install (TARGETS ActsJsonPlugin
         EXPORT ActsJsonPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsJsonPlugin TARGETS ActsJsonPlugin)
