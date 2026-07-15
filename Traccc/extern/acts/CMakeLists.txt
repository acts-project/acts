# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.25 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Acts (including Detray) as part of the TRACCC project" )

# Declare where to get Acts from.
set( TRACCC_ACTS_SOURCE
   GIT_REPOSITORY https://github.com/acts-project/acts.git
   GIT_TAG cda8b940cedfaa413931184587d6ac200e6198d0
   CACHE STRING "Source for Acts, when built as part of this project" )
mark_as_advanced( TRACCC_ACTS_SOURCE )
FetchContent_Declare( Acts SYSTEM ${TRACCC_ACTS_SOURCE} )

# Options used in the build of Acts.
set( ACTS_SETUP_EIGEN3 TRUE CACHE BOOL
   "Set up Eigen in the Acts code" )
set( ACTS_USE_SYSTEM_EIGEN3 FALSE CACHE BOOL
   "Build Eigen3 from source" )
set( ACTS_BUILD_PLUGIN_JSON TRUE CACHE BOOL
   "Build JSON plugin in Acts" )
set( ACTS_SETUP_DETRAY TRUE CACHE BOOL
   "Set up Detray in the Acts code" )

# Set the detray log level
set( DETRAY_SET_LOGGING "${TRACCC_DEVICE_LOG_LVL}" CACHE STRING
   "Log level to use in device code" FORCE)
message( STATUS "Setting DETRAY and TRACCC device log level to ${DETRAY_SET_LOGGING}" )

# Test utils and validation tools are needed for traccc unittests/ACTS detray plugin
set(DETRAY_BUILD_VALIDATION_TOOLS TRUE CACHE BOOL "")
set(DETRAY_BUILD_TEST_UTILS TRUE CACHE BOOL "")
# HACK: Make sure that ACTS builds ACTSVG, which is required by Detray validation tools.
set( ACTS_BUILD_PLUGIN_ACTSVG TRUE CACHE BOOL "" )

# The ITk metadata is required, so explicitly request that it be built.
if (NOT("itk_metadata" IN_LIST DETRAY_GENERATE_METADATA))
   list(APPEND DETRAY_GENERATE_METADATA "itk_metadata")
endif()

# Get it into the current directory.
FetchContent_MakeAvailable( Acts )
