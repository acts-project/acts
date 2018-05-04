file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsMaterialPlugin SHARED ${src_files})

target_include_directories (ActsMaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ActsMaterialPlugin PUBLIC ActsCore ActsLegacy)

install (TARGETS ActsMaterialPlugin
         EXPORT ActsMaterialPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT MaterialPlugin TARGETS ActsMaterialPlugin)
