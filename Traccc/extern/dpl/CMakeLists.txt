# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.25 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building oneDPL as part of the TRACCC project" )

# Declare where to get DPL from.
set( TRACCC_DPL_SOURCE
   "URL;https://github.com/uxlfoundation/oneDPL/archive/refs/tags/oneDPL-2022.11.0-release.tar.gz;URL_MD5;d42f971a7fb53574e511e207b22151cd;PATCH_COMMAND;patch;-p1;-i;${CMAKE_CURRENT_SOURCE_DIR}/oneDPL-2022.11.0.patch"
   CACHE STRING "Source for DPL, when built as part of this project" )
mark_as_advanced( TRACCC_DPL_SOURCE )
# Must not use SYSTEM here, as the oneAPI compiler then ignores this version
# of oneDPL.
FetchContent_Declare( DPL ${TRACCC_DPL_SOURCE} )

# Set the default oneDPL threading backend.
set( ONEDPL_BACKEND "dpcpp" CACHE STRING "oneDPL threading backend" )
set( SYCL_SUPPORT ON )

# Get it into the current directory.
FetchContent_MakeAvailable( DPL )
