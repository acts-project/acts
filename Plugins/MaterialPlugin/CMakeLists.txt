file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (MaterialPlugin SHARED ${src_files})

target_include_directories (MaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (MaterialPlugin PUBLIC Acts::Core Acts::Legacy)

install (TARGETS MaterialPlugin
         EXPORT ActsMaterialPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT MaterialPlugin TARGETS MaterialPlugin)
