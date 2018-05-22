file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ACTSDD4hepPlugin SHARED ${src_files})

target_include_directories (ACTSDD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories (ACTSDD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})
target_link_libraries (ACTSDD4hepPlugin PUBLIC Acts::Core)
target_link_libraries (ACTSDD4hepPlugin PUBLIC ACTSTGeoPlugin)
target_link_libraries (ACTSDD4hepPlugin PUBLIC ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDSEGMENTATION_LIBRARY})

install (TARGETS ACTSDD4hepPlugin
         EXPORT ACTSDD4hepPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT DD4hepPlugin TARGETS ACTSDD4hepPlugin)
