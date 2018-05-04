file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library (ActsTGeoPlugin SHARED ${src_files})

target_include_directories(ActsTGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsTGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ActsTGeoPlugin PUBLIC ActsCore)
target_link_libraries(ActsTGeoPlugin PUBLIC ROOT::Geom)

install (TARGETS ActsTGeoPlugin EXPORT ActsTGeoPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT TGeoPlugin TARGETS ActsTGeoPlugin)
