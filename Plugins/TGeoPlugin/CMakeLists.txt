find_package(ROOT 6.10 REQUIRED COMPONENTS Geom)

file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ACTSTGeoPlugin SHARED ${src_files})

target_include_directories(ACTSTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTSTGeoPlugin PUBLIC ACTSCore)
target_link_libraries(ACTSTGeoPlugin PUBLIC ROOT::Geom)

install (TARGETS ACTSTGeoPlugin EXPORT ACTSTGeoPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT TGeoPlugin TARGETS ACTSTGeoPlugin)
