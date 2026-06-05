# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.25 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the TRACCC project" )

# Declare where to get dfelibs from.
set( TRACCC_DFELIBS_SOURCE
   "URL;https://github.com/acts-project/dfelibs/archive/refs/tags/v20211029.tar.gz;URL_MD5;87fb09c5a11b98250f5e266e9cd501ea"
   CACHE STRING "Source for dfelibs, when built as part of this project" )
mark_as_advanced( TRACCC_DFELIBS_SOURCE )
FetchContent_Declare( dfelibs SYSTEM ${TRACCC_DFELIBS_SOURCE} )

# Options used in the build of dfelibs.
set( dfelibs_BUILD_EXAMPLES FALSE CACHE BOOL
   "Turn off the build of the dfelibs examples" )
set( dfelibs_BUILD_UNITTESTS FALSE CACHE BOOL
   "Turn off the build of the dfelibs unit tests" )
set( dfelibs_ENABLE_INSTALL TRUE CACHE BOOL
   "Have dfelibs be installed together with the rest of the project" )

# Make dfelibs work without warnings with modern CMake versions.
set( CMAKE_POLICY_VERSION_MINIMUM "3.10" )

# Get it into the current directory.
FetchContent_MakeAvailable( dfelibs )
