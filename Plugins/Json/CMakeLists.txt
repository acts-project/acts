file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

acts_add_library (JsonPlugin SHARED ${src_files})

target_include_directories (JsonPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries (JsonPlugin PUBLIC Acts::Core)

install (TARGETS JsonPlugin
         EXPORT ActsJsonTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

acts_add_targets_to_cdash_project(PROJECT JsonPlugin TARGETS JsonPlugin)
