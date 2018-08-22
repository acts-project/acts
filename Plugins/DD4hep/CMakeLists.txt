file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (DD4hepPlugin SHARED ${src_files})

target_include_directories (DD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories (DD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})
target_link_libraries (DD4hepPlugin PUBLIC Acts::Core)
target_link_libraries (DD4hepPlugin PUBLIC Acts::TGeoPlugin Acts::DigitizationPlugin)
target_link_libraries (DD4hepPlugin PUBLIC ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDSEGMENTATION_LIBRARY})

install (TARGETS DD4hepPlugin
         EXPORT ActsDD4hepTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT DD4hepPlugin TARGETS DD4hepPlugin)
