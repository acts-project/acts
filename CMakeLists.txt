cmake_minimum_required (VERSION 3.5)

project (ACTS VERSION 0.6.0 LANGUAGES CXX)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# build options

# all options must be defined here to keep the list directly visible and
# make the option variables available everywhere
option (BUILD_DOC "Build documentation" OFF)
option (BUILD_EXAMPLES "Build examples" ON)
option (BUILD_TESTS "Build unit tests" ON)
option (BUILD_INTEGRATION_TESTS "Build integration tests" OFF)
option (BUILD_DD4HEP_PLUGIN "Build DD4HEP plugins" OFF)
option (BUILD_MATERIAL_PLUGIN "Build Material plugins" ON)
option (BUILD_TGEO_PLUGIN "Build TGeo plugins" OFF)
option (SUPPRESS_LIBRARY_WARNINGS "Suppress warnings from Eigen and Boost libraries" OFF)
# all other compile-time parameters must be defined here for clear visibility
# and to avoid forgotten options somewhere deep in the hierarchy
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "ACTS/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "Default track parameter definition")

# handle inter-plugin dependencies
if (BUILD_DD4HEP_PLUGIN)
  set (BUILD_TGEO_PLUGIN ON)
endif()

# required packages

find_package (Boost 1.62 REQUIRED)
find_package (Eigen 3.2.9 REQUIRED)

# additional configuration and tools

include (GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include (ACTSCompilerOptions) # default compile options
include (ACTSFunctions) # ACTS-specific helper functions
include (ACTSComponentsHelpers) # handle components via add_..._if commands

# Configure CCache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  message(STATUS "Use ccache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# silence warning about missing RPATH on Mac OSX
set (CMAKE_MACOSX_RPATH 1)

# convert the version into a plain number to be used in the version header file
string (REPLACE "." "" ACTS_VERSION_NUMBER ${PROJECT_VERSION})

# core library and plugins
add_component (Core)
add_subdirectory (Plugins)

# examples and tests
if (BUILD_TESTS OR BUILD_INTEGRATION_TESTS)
  include (CTest) # must be included from the main CMakeLists.txt
endif ()
add_component (Examples)
add_subdirectory_if (Tests BUILD_TESTS)
add_subdirectory_if (IntegrationTests BUILD_INTEGRATION_TESTS)

# documentation
add_subdirectory_if (doc BUILD_DOC)

# add version header
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/ACTSVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTSVersion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ACTS)

# add setup scripts
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION ${CMAKE_INSTALL_BINDIR})

# add CMake infrastructure
include (ACTSCreateConfig)
