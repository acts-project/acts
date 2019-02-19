cmake_minimum_required(VERSION 3.5)

# read version number from file
file(READ version_number _acts_version)
string(STRIP ${_acts_version} _acts_version)

project (Acts VERSION ${_acts_version} LANGUAGES CXX)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# get git commit hash
include (GetGitRevisionDescription)
get_git_head_revision(_git_refspec _git_hash)
git_local_changes(_git_local_changes)
string(SUBSTRING "${_git_hash}" 0 8 _git_hash_short)

# build options

# all options must be defined here to keep the list directly visible and
# make the option variables available everywhere
option (ACTS_BUILD_LEGACY "Build legacy package" OFF)
option (ACTS_BUILD_DOC "Build documentation" OFF)
option (ACTS_BUILD_EXAMPLES "Build examples" OFF)
option (ACTS_BUILD_TESTS "Build unit tests" OFF)
option (ACTS_BUILD_INTEGRATION_TESTS "Build integration tests" OFF)
option (ACTS_BUILD_DIGITIZATION_PLUGIN "Build Digitization plugin" OFF)
option (ACTS_BUILD_DD4HEP_PLUGIN "Build DD4HEP plugin" OFF)
option (ACTS_BUILD_IDENTIFICATION_PLUGIN "Build Identification plugin" OFF)
option (ACTS_BUILD_JSON_PLUGIN "Build Json plugin" OFF)
option (ACTS_BUILD_MATERIAL_PLUGIN "Build Material plugin" OFF)
option (ACTS_BUILD_TGEO_PLUGIN "Build TGeo plugin" OFF)

# all other compile-time parameters must be defined here for clear visibility
# and to avoid forgotten options somewhere deep in the hierarchy
set (ACTS_PARAMETER_DEFINITIONS_PLUGIN "Acts/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "Default track parameter definition")

# handle inter-plugin dependencies
# DD4hepPlugin depends on TGeoPlugin
if (ACTS_BUILD_DD4HEP_PLUGIN)
  set (ACTS_BUILD_TGEO_PLUGIN ON)
endif()
# TGeoPlugin depends on IdentificationPlugin
if (ACTS_BUILD_TGEO_PLUGIN)
  set (ACTS_BUILD_IDENTIFICATION_PLUGIN ON)
endif()

# additional configuration and tools
include (GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include (ActsCompilerOptions) # default compile options
include (ActsFunctions) # Acts-specific helper functions
include (ActsComponentsHelpers) # handle components via add_..._if commands

# required packages
find_package(Boost 1.62 REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Eigen 3.2.9 REQUIRED)

# optional packages
if(ACTS_BUILD_DOC)
  find_package(Doxygen 1.8.11 REQUIRED)
endif()
if(ACTS_BUILD_DD4HEP_PLUGIN)
  # Work around DD4hep forcing CMAKE_CXX_STANDARD to 14
  set(_cmake_cxx_standard ${CMAKE_CXX_STANDARD})
  find_package(DD4hep 1.2 REQUIRED COMPONENTS DDCore)
  set(CMAKE_CXX_STANDARD ${_cmake_cxx_standard})
endif()
if(ACTS_BUILD_TGEO_PLUGIN)
  find_package(ROOT 6.10 REQUIRED COMPONENTS Geom)
endif()

# core library and plugins
add_component (Core Core)
add_component_if (Legacy Legacy ACTS_BUILD_LEGACY)
add_subdirectory (Plugins)

# examples and tests
if(ACTS_BUILD_TESTS OR ACTS_BUILD_INTEGRATION_TESTS)
  include(CTest) # must be included from the main CMakeLists.txt
endif()
add_subdirectory_if (Tests ACTS_BUILD_TESTS OR ACTS_BUILD_INTEGRATION_TESTS OR ACTS_BUILD_EXAMPLES)

# documentation
add_subdirectory_if (doc ACTS_BUILD_DOC)

# setup script and cmake configuration files
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION ${CMAKE_INSTALL_BINDIR})
include (ActsCreateConfig)
