# get list of all source files
file (GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

# define library target
add_library (ACTSMaterialPlugin SHARED ${src_files})

# setup include directories
target_include_directories (ACTSMaterialPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

# setup linker dependencies
target_link_libraries (ACTSMaterialPlugin PUBLIC ACTSCore)

# set installation directories
install (TARGETS ACTSMaterialPlugin
         EXPORT ACTSMaterialPluginTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY include/ACTS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# add to CDash MaterialPlugin subproject
acts_add_targets_to_cdash_project(PROJECT MaterialPlugin TARGETS ACTSMaterialPlugin)
