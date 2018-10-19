file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsDD4hepPlugin SHARED ${src_files})

target_include_directories (ActsDD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories (ActsDD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})
target_link_libraries (ActsDD4hepPlugin PUBLIC ActsCore)
target_link_libraries (ActsDD4hepPlugin PUBLIC ActsTGeoPlugin)
target_link_libraries (ActsDD4hepPlugin PUBLIC ${DD4hep_DDCORE_LIBRARY})

install (TARGETS ActsDD4hepPlugin
         EXPORT ActsDD4hepPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT ActsDD4hepPlugin TARGETS ActsDD4hepPlugin)
