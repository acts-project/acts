# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.25 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Fetching covfie as part of the traccc project" )

# Declare where to get covfie from.
set( TRACCC_COVFIE_SOURCE
   "URL;https://github.com/acts-project/covfie/archive/refs/tags/v0.15.4.tar.gz;URL_MD5;80de9b43644c59d8239fb8cc184c6276"
   CACHE STRING "Source for covfie, when built as part of this project" )
mark_as_advanced( TRACCC_COVFIE_SOURCE )
FetchContent_Declare( covfie SYSTEM ${TRACCC_COVFIE_SOURCE} )

# Options used for covfie.
set( COVFIE_BUILD_EXAMPLES OFF CACHE BOOL "Build covfie examples")
set( COVFIE_BUILD_TESTS OFF CACHE BOOL "Build covfie tests")
set( COVFIE_BUILD_BENCHMARKS OFF CACHE BOOL "Build covfie benchmarks")

set( COVFIE_PLATFORM_CPU ON CACHE BOOL "Enable covfie CPU platform")
set( COVFIE_PLATFORM_CUDA ${TRACCC_BUILD_CUDA} CACHE BOOL "Enable covfie CUDA platform")
set( COVFIE_PLATFORM_SYCL ${TRACCC_BUILD_SYCL} CACHE BOOL "Enable covfie SYCL platform")

set( COVFIE_QUIET ON CACHE BOOL "Quiet covfie feature warnings")

# Get it into the current directory.
FetchContent_MakeAvailable( covfie )
