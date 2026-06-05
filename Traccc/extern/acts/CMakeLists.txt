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
FetchContent_Declare( Acts SYSTEM )

# Options used in the build of Acts.
set( ACTS_SETUP_EIGEN3 TRUE CACHE BOOL
   "Set up Eigen in the Acts code" )
set( ACTS_USE_SYSTEM_EIGEN3 FALSE CACHE BOOL
   "Build Eigen3 from source" )
set( ACTS_BUILD_PLUGIN_JSON TRUE CACHE BOOL
   "Build JSON plugin in Acts" )
set( ACTS_SETUP_DETRAY TRUE CACHE BOOL
   "Set up Detray in the Acts code" )
set(DETRAY_BUILD_TEST_UTILS TRUE CACHE BOOL "")
# The ITk metadata is required, so explicitly request that it be built.
if (NOT("itk_metadata" IN_LIST DETRAY_GENERATE_METADATA))
   list(APPEND DETRAY_GENERATE_METADATA "itk_metadata")
endif()

FetchContent_MakeAvailable( Acts )
