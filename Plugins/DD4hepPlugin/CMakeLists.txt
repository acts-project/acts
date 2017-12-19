find_package (DD4hep ${ACTS_FIND_QUIET_TEXT} REQUIRED COMPONENTS DDCore DDSegmentation)

# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

# define library target
add_library (ACTSDD4hepPlugin SHARED ${src_files})

# setup include directories
target_include_directories (ACTSDD4hepPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories (ACTSDD4hepPlugin PUBLIC ${DD4hep_INCLUDE_DIRS})

# setup linker dependencies
target_link_libraries (ACTSDD4hepPlugin PUBLIC ACTSCore)
target_link_libraries (ACTSDD4hepPlugin PUBLIC ACTSTGeoPlugin)
target_link_libraries (ACTSDD4hepPlugin PUBLIC ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDSEGMENTATION_LIBRARY})

# set installation directories
install (TARGETS ACTSDD4hepPlugin
         EXPORT ACTSDD4hepPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# add to CDash DD4hepPlugin subproject
acts_add_targets_to_cdash_project(PROJECT DD4hepPlugin TARGETS ACTSDD4hepPlugin)
