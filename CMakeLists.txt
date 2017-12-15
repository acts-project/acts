cmake_minimum_required (VERSION 3.5)

project (ACTS VERSION 0.6.0 LANGUAGES CXX)

# all options must be defined here to keep the list directly visible and
# make the option variables available everywhere
option (BUILD_DOC "Build documentation" OFF)
option (BUILD_TESTS "Build unit tests" ON)
option (BUILD_INTEGRATION_TESTS "Build integration tests" OFF)
option (BUILD_DD4HEP_PLUGIN "build DD4HEP plugins" OFF)
option (BUILD_MATERIAL_PLUGIN "build Material plugins" ON)
option (BUILD_TGEO_PLUGIN "build TGeo plugins" OFF)
option (SUPPRESS_LIBRARY_WARNINGS "Suppress warnings from Eigen and Boost libraries" OFF)

# handle inter-plugin dependencies
if (BUILD_DD4HEP_PLUGIN)
  set (BUILD_TGEO_PLUGIN ON)
endif()

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# print all available options
if (DEFINED PRINT_OPTIONS)
  # silence cmake
  set (_print_options ${PRINT_OPTIONS})
  # include options file
  include (cmake/PrintOptions.cmake)
  unset (PRINT_OPTIONS CACHE)
  return ()
endif (DEFINED PRINT_OPTIONS)

# include support for colored output
include (cmake/Colors.cmake)

# Use GNU-style hierarchy for installing build products
include (GNUInstallDirs)

# include ACTS specific compilation options
include (cmake/ACTSCompilerOptions.cmake)

# include ACTS specific macros
include (cmake/ACTSFunctions.cmake)

# include versions for dependencies of ACTS
include (cmake/ACTSDependencies.cmake)

# Configure CCache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message (STATUS "${Yellow}use ccache${ColorReset}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# silence warning about missing RPATH on Mac OSX
set (CMAKE_MACOSX_RPATH 1)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${PROJECT_VERSION})

# default parameter definitions
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "default track parameter definition")

# finding boost on lxplus does not work with standard boost-cmake look up
set (Boost_NO_BOOST_CMAKE ON CACHE BOOL "disable boost-cmake find mechanism")

add_subdirectory (Core)
add_subdirectory (doc)
add_subdirectory (Plugins)
add_subdirectory (Examples)

# add unit tests
include (CTest)
add_subdirectory (Tests)
add_subdirectory (IntegrationTests)

# add version header
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add setup scripts
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION ${CMAKE_INSTALL_BINDIR})

# list components which are going to be installed
message (STATUS "going to install the following components:")
foreach (_scomp ${_supported_components})
  message (STATUS "   - ${_scomp}")
endforeach ()

# add CMake infrastructure
include (cmake/ACTSCreateConfig.cmake)
