find_package(ROOT 6.10 REQUIRED COMPONENTS Geom)

file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ACTSTGeoPlugin SHARED ${src_files})

target_include_directories (ACTSTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories (ACTSTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (ACTSTGeoPlugin PUBLIC ACTSCore)
# temporary support for ROOT version 6.08 without CMake Namespace
if (${ROOT_VERSION_MINOR} MATCHES "8" )
  target_link_libraries (ACTSTGeoPlugin PUBLIC Geom)
else (${ROOT_VERSION_MINOR} MATCHES "8" )
  target_link_libraries (ACTSTGeoPlugin PUBLIC ROOT::Geom)
endif (${ROOT_VERSION_MINOR} MATCHES "8" )

install (TARGETS ACTSTGeoPlugin EXPORT ACTSTGeoPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT TGeoPlugin TARGETS ACTSTGeoPlugin)
