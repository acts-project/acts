cmake_minimum_required(VERSION 3.11)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# determine project version; sets _acts_version and _acts_commit_hash
include(ActsRetrieveVersion)

# Beomki: Need to modulize the CUDA activation later....
project(Acts VERSION ${_acts_version} LANGUAGES CUDA CXX)			  
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_FLAGS -arch=sm_61)


# build options

# all options must be defined here to keep the list directly visible and
# make the option variables available everywhere
option(ACTS_BUILD_DD4HEP_PLUGIN "Build DD4HEP plugin" OFF)
option(ACTS_BUILD_DIGITIZATION_PLUGIN "Build Digitization plugin" OFF)
option(ACTS_BUILD_IDENTIFICATION_PLUGIN "Build Identification plugin" OFF)
option(ACTS_BUILD_JSON_PLUGIN "Build Json plugin" OFF)
option(ACTS_USE_BUNDLED_NLOHMANN_JSON "Use external nlohmann::json" OFF)
option(ACTS_BUILD_TGEO_PLUGIN "Build TGeo plugin" OFF)
option(ACTS_BUILD_FATRAS "Build FAst TRAcking Simulation package" OFF)
option(ACTS_BUILD_LEGACY "Build legacy package" OFF)
option(ACTS_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(ACTS_BUILD_EXAMPLES "Build examples" OFF)
option(ACTS_BUILD_UNITTESTS "Build unit tests" ON)
option(ACTS_BUILD_INTEGRATIONTESTS "Build integration tests" OFF)
option(ACTS_BUILD_DOC "Build documentation" OFF)
# all other compile-time parameters must be defined here for clear visibility
# and to avoid forgotten options somewhere deep in the hierarchy
set(ACTS_PARAMETER_DEFINITIONS_HEADER "" CACHE FILEPATH "Use a different (track) parameter definitions header")

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
find_package(Boost 1.69 REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Eigen 3.2.9 REQUIRED)

# optional packages
if(ACTS_BUILD_DD4HEP_PLUGIN)
  find_package(DD4hep 1.2 REQUIRED COMPONENTS DDCore)
endif()
if(ACTS_BUILD_JSON_PLUGIN)
  if(ACTS_USE_BUNDLED_NLOHMANN_JSON)
    message(STATUS "Will build nlohmann::json as part of this build")
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(thirdparty/nlohmann_json)
  else()
    message(STATUS "Using external install of nlohmann::json")
    find_package(nlohmann_json 3.2.0 REQUIRED)
  endif()
endif()
if(ACTS_BUILD_TGEO_PLUGIN)
  find_package(ROOT 6.10 REQUIRED COMPONENTS Geom)
endif()
if(ACTS_BUILD_DOC)
  find_package(Doxygen 1.8.11 REQUIRED)
endif()

# core library, core plugins, and other components
add_component(Core Core)
add_subdirectory(Plugins)
add_component(Fatras Fatras)

# benchmarks, examples, and tests
if(ACTS_BUILD_UNITTESTS OR ACTS_BUILD_INTEGRATIONTESTS)
  enable_testing() # must be set in the main CMakeLists.txt
endif()
add_subdirectory(Tests)

# documentation
add_subdirectory_if(doc ACTS_BUILD_DOC)

# create cmake configuration files and environment setup script
include(ActsCreateCMakeConfig)
include(ActsCreateSetup)
