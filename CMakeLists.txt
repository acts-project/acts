cmake_minimum_required(VERSION 3.11)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# determine project version; sets _acts_version and _acts_commit_hash
include(ActsRetrieveVersion)

project(Acts VERSION ${_acts_version} LANGUAGES CXX)

# build options

# all options must be defined here to keep the list directly visible and
# make the option variables available everywhere
option(ACTS_BUILD_LEGACY "Build legacy package" OFF)
option(ACTS_BUILD_DOC "Build documentation" OFF)
option(ACTS_BUILD_EXAMPLES "Build examples" OFF)
option(ACTS_BUILD_TESTS "Build unit tests" OFF)
option(ACTS_BUILD_INTEGRATION_TESTS "Build integration tests" OFF)
option(ACTS_BUILD_DD4HEP_PLUGIN "Build DD4HEP plugin" OFF)
option(ACTS_BUILD_DIGITIZATION_PLUGIN "Build Digitization plugin" OFF)
option(ACTS_BUILD_IDENTIFICATION_PLUGIN "Build Identification plugin" OFF)
option(ACTS_BUILD_JSON_PLUGIN "Build Json plugin" OFF)
option(ACTS_USE_EXTERNAL_NLOHMANN_JSON "Use external nlohmann::json" OFF)
option(ACTS_BUILD_TGEO_PLUGIN "Build TGeo plugin" OFF)

# all other compile-time parameters must be defined here for clear visibility
# and to avoid forgotten options somewhere deep in the hierarchy
set(ACTS_PARAMETER_DEFINITIONS_PLUGIN "Acts/Utilities/detail/DefaultParameterDefinitions.hpp" CACHE FILEPATH "Default track parameter definition")

# handle inter-plugin dependencies
# DD4hepPlugin depends on TGeoPlugin
if(ACTS_BUILD_DD4HEP_PLUGIN)
  set(ACTS_BUILD_TGEO_PLUGIN ON)
endif()
# TGeoPlugin depends on IdentificationPlugin
if(ACTS_BUILD_TGEO_PLUGIN)
  set(ACTS_BUILD_IDENTIFICATION_PLUGIN ON)
endif()

# additional configuration and tools
include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include(ActsCompilerOptions) # default compile options
include(ActsComponentsHelpers) # handle components via add_..._if commands
include(TargetSourcesLocal) # needed only until we require CMake 3.13

# required packages
set(Boost_NO_BOOST_CMAKE ON) # disable new cmake features from Boost 1.70 on
find_package(Boost 1.62 REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Eigen 3.2.9 REQUIRED)

# optional packages
if(ACTS_BUILD_DOC)
  find_package(Doxygen 1.8.11 REQUIRED)
endif()
if(ACTS_BUILD_DD4HEP_PLUGIN)
  find_package(DD4hep 1.2 REQUIRED COMPONENTS DDCore)
endif()
if(ACTS_BUILD_TGEO_PLUGIN)
  find_package(ROOT 6.10 REQUIRED COMPONENTS Geom)
endif()

# core library and plugins
add_component(Core Core)
add_component_if(Legacy Legacy ACTS_BUILD_LEGACY)
add_subdirectory(Plugins)
# examples and tests
if(ACTS_BUILD_TESTS OR ACTS_BUILD_INTEGRATION_TESTS)
  enable_testing() # must be set from the main CMakeLists.txt
endif()
add_subdirectory_if(Tests ACTS_BUILD_TESTS OR ACTS_BUILD_INTEGRATION_TESTS OR ACTS_BUILD_EXAMPLES)
# documentation
add_subdirectory_if(doc ACTS_BUILD_DOC)

# create cmake configuration files and environment setup script
include(ActsCreateCMakeConfig)
include(ActsCreateSetup)
