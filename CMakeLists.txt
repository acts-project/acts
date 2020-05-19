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

# all options and compile-time parameters must be defined here for clear
# visibility and make them available everywhere
# core related options
set(ACTS_PARAMETER_DEFINITIONS_HEADER "" CACHE FILEPATH "Use a different (track) parameter definitions header")
option(ACTS_BUILD_DD4HEP_PLUGIN "Build DD4HEP plugin" OFF)
option(ACTS_BUILD_DIGITIZATION_PLUGIN "Build Digitization plugin" OFF)
option(ACTS_BUILD_IDENTIFICATION_PLUGIN "Build Identification plugin" OFF)
option(ACTS_BUILD_JSON_PLUGIN "Build Json plugin" OFF)
option(ACTS_USE_BUNDLED_NLOHMANN_JSON "Use external nlohmann::json" ON)
option(ACTS_BUILD_TGEO_PLUGIN "Build TGeo plugin" OFF)
option(ACTS_BUILD_CUDA_PLUGIN "Build Cuda plugin" OFF)
option(ACTS_BUILD_LEGACY "Build legacy plugin" OFF)
option(ACTS_BUILD_FATRAS "Build FAst TRAcking Simulation package" OFF)
# examples related options
option(ACTS_BUILD_EXAMPLES "Build standalone examples" OFF)
option(ACTS_BUILD_EXAMPLES_DD4HEP "Build DD4hep-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_GEANT4 "Build Geant4-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_HEPMC3 "Build HepMC3-based code in the examples" OFF)
option(ACTS_BUILD_EXAMPLES_PYTHIA8 "Build Pythia8-based code in the examples" OFF)
# test related options
option(ACTS_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(ACTS_BUILD_UNITTESTS "Build unit tests" OFF)
option(ACTS_BUILD_INTEGRATIONTESTS "Build integration tests" OFF)
# other options
option(ACTS_BUILD_EVERYTHING "Build with most options enabled (except documentation)" OFF)
option(ACTS_BUILD_DOCS "Build documentation" OFF)

# handle option inter-dependencies and the everything flag
# DD4hep and HepMC3 are excluded from the EVERYTHING option for now since
# they are not yet available in every CI build.
include(ActsOptionHelpers)
# options without inter-dependencies
set_option_if(ACTS_BUILD_BENCHMARKS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_UNITTESTS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_INTEGRATIONTESTS ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_DD4HEP ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_GEANT4 ACTS_BUILD_EVERYTHING)
# set_option_if(ACTS_BUILD_EXAMPLES_HEPMC3 ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_EXAMPLES_PYTHIA8 ACTS_BUILD_EVERYTHING)
# any examples component activates the general option
set_option_if(
  ACTS_BUILD_EXAMPLES
  ACTS_BUILD_EXAMPLES_DD4HEP
  OR ACTS_BUILD_EXAMPLES_GEANT4
  OR ACTS_BUILD_EXAMPLES_HEPMC3
  OR ACTS_BUILD_EXAMPLES_PYTHIA8
  OR ACTS_BUILD_EVERYTHING)
# core plugins might be required by examples or depend on each other
set_option_if(
  ACTS_BUILD_DD4HEP_PLUGIN
  ACTS_BUILD_EXAMPLES_DD4HEP) # OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_TGEO_PLUGIN
  ACTS_BUILD_DD4HEP_PLUGIN OR ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_IDENTIFICATION_PLUGIN
  ACTS_BUILD_TGEO_PLUGIN OR ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_DIGITIZATION_PLUGIN
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_JSON_PLUGIN
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(
  ACTS_BUILD_FATRAS
  ACTS_BUILD_EXAMPLES OR ACTS_BUILD_EVERYTHING)
set_option_if(ACTS_BUILD_LEGACY ACTS_BUILD_EVERYTHING)

# additional configuration and tools
include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include(ActsCompilerOptions) # default compile options
include(ActsComponentsHelpers) # handle components via add_..._if commands
include(TargetSourcesLocal) # needed only until we require CMake 3.13
# place build products in `<build_dir>/bin` and `<build_dir>/lib` for simple use
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# required packages
set(Boost_NO_BOOST_CMAKE ON) # disable new cmake features from Boost 1.70 on
find_package(Boost 1.69 REQUIRED COMPONENTS filesystem program_options unit_test_framework)
find_package(Eigen 3.2.9 REQUIRED)

# optional packages
#
# find packages explicitely for each component even if this can mean
# searching for the same package twice. This avoid having complex if/else
# trees to sort of when a particular package is actually needed.
#
# version variables for external dependencies that are required more than once
set(acts_dd4hep_version 1.10)
set(acts_root_version 6.10)
if(ACTS_BUILD_DD4HEP_PLUGIN)
  find_package(DD4hep ${acts_dd4hep_version} REQUIRED COMPONENTS DDCore)
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
  find_package(ROOT ${acts_root_version} REQUIRED COMPONENTS Geom)
endif()
if(ACTS_BUILD_EXAMPLES)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  # we could select ROOT components based on the configured core plugins and
  # standalone components. for simplicity always request all possible
  # required components.
  find_package(ROOT ${acts_root_version} REQUIRED COMPONENTS Core Geom GenVector Hist Tree TreePlayer)
  find_package(TBB REQUIRED)
  add_subdirectory(thirdparty/dfelibs)
endif()
if(ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_GEANT4)
  find_package(DD4hep ${acts_dd4hep_version} REQUIRED COMPONENTS DDCore DDG4)
elseif(ACTS_BUILD_EXAMPLES_DD4HEP)
  find_package(DD4hep ${acts_dd4hep_version} REQUIRED COMPONENTS DDCore)
endif()
if(ACTS_BUILD_EXAMPLES_GEANT4)
  find_package(Geant4 REQUIRED COMPONENTS gdml)
endif()
if(ACTS_BUILD_EXAMPLES_HEPMC3)
  # must be version 3.1+. version 3.0 uses HepMC w/o number as module name.
  # version 3.1.0 misses the version variable so we can not check for it
  find_package(HepMC3 REQUIRED)
  find_package(HepPID REQUIRED)
  # HepMC version 3.1.0 has a broken cmake config: it searches for the
  # HepMC libray instead of the correct HepMC3 one and sets HEPMC_LIBRARIES
  # instead of the expected HEPMC3_LIBRARIES. fix manually.
  if(NOT HEPMC3_LIBRARIES)
    message(STATUS "Workaround for HepMC3. Search in $ENV{HepMC3_DIR}")
    find_library(
      HEPMC3_LIB Names HepMC3
      PATHS $ENV{HepMC3_DIR}/lib $ENV{HepMC3_DIR}/lib64 NO_DEFAULT_PATH)
    set(HEPMC3_INCLUDE_DIR ${HEPMC_INCLUDE_DIR})
    set(HEPMC3_LIBRARIES ${HEPMC3_LIB})
  endif()
endif()
if(ACTS_BUILD_EXAMPLES_PYTHIA8)
  find_package(Pythia8 REQUIRED)
endif()
if(ACTS_BUILD_DOCS)
  find_package(Doxygen 1.8.11 REQUIRED)
  find_package(Sphinx REQUIRED)
endif()
if(ACTS_BUILD_CUDA_PLUGIN)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# core library, core plugins, and other components
add_component(Core Core)
add_subdirectory(Plugins)
add_component_if(Fatras Fatras ACTS_BUILD_FATRAS)

# manual examples
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)

# automated tests and benchmarks
if(ACTS_BUILD_UNITTESTS OR ACTS_BUILD_INTEGRATIONTESTS)
  enable_testing() # must be set in the main CMakeLists.txt
endif()
add_subdirectory(Tests)

# documentation
add_subdirectory_if(docs ACTS_BUILD_DOCS)

# create cmake configuration files and environment setup script
include(ActsCreateCMakeConfig)
include(ActsCreateSetup)
