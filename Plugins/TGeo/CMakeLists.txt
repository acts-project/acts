file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (TGeoPlugin SHARED ${src_files})

target_include_directories(TGeoPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(TGeoPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(TGeoPlugin PUBLIC Acts::Core Acts::DigitizationPlugin)
target_link_libraries(TGeoPlugin PUBLIC ROOT::Geom)

install (TARGETS TGeoPlugin EXPORT ActsTGeoTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT TGeoPlugin TARGETS TGeoPlugin)
