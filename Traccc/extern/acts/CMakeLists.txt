# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.25)

# Tell the user what's happening.
message(STATUS "Building Acts (including Detray) as part of the TRACCC project")

# Acts is colocated in this repository: Traccc lives in-tree under the Acts
# source, which is the single source of truth. Build against that copy directly
# (the repository root, three levels up from here) rather than fetching an
# upstream copy. The add_subdirectory() happens at the end of this file, once the
# Acts/Detray build options below have been set.
get_filename_component(
    TRACCC_ACTS_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../../.."
    ABSOLUTE
)

# Options used in the build of Acts.
set(ACTS_SETUP_EIGEN3 TRUE CACHE BOOL "Set up Eigen in the Acts code")
set(ACTS_USE_SYSTEM_EIGEN3 FALSE CACHE BOOL "Build Eigen3 from source")
set(ACTS_BUILD_PLUGIN_JSON TRUE CACHE BOOL "Build JSON plugin in Acts")
set(ACTS_SETUP_DETRAY TRUE CACHE BOOL "Set up Detray in the Acts code")

# Acts and Detray are in-tree, first-party code here, so build them with
# warnings-as-errors (matching Acts' own CI). This catches warnings in their own
# compilation -- not just those in headers pulled into traccc translation units,
# which are already covered by dropping the -isystem treatment.
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

# Set the detray log level
set(DETRAY_SET_LOGGING
    "${TRACCC_DEVICE_LOG_LVL}"
    CACHE STRING
    "Log level to use in device code"
    FORCE
)
message(
    STATUS
    "Setting DETRAY and TRACCC device log level to ${DETRAY_SET_LOGGING}"
)

# Test utils and validation tools are needed for traccc unittests/ACTS detray plugin
set(DETRAY_BUILD_VALIDATION_TOOLS TRUE CACHE BOOL "")
set(DETRAY_BUILD_TEST_UTILS TRUE CACHE BOOL "")
# HACK: Make sure that ACTS builds ACTSVG, which is required by Detray validation tools.
set(ACTS_BUILD_PLUGIN_ACTSVG TRUE CACHE BOOL "")

# The ITk metadata is required, so explicitly request that it be built.
if(NOT ("itk_metadata" IN_LIST DETRAY_GENERATE_METADATA))
    list(APPEND DETRAY_GENERATE_METADATA "itk_metadata")
endif()

# Build the colocated Acts (and, transitively, Detray) as a subdirectory of this
# build. An explicit binary dir is required because the source lives outside the
# current directory tree.
add_subdirectory(
    "${TRACCC_ACTS_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/acts"
)
