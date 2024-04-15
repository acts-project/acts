# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building VecMem as part of the Acts project" )

set( VECMEM_VERSION "v${_acts_vecmem_version}")

# Declare where to get VecMem from.
set( ACTS_VECMEM_GIT_REPOSITORY "https://github.com/acts-project/vecmem.git"
   CACHE STRING "Git repository to take VecMem from" )
set( ACTS_VECMEM_GIT_TAG "${VECMEM_VERSION}" CACHE STRING "Version of VecMem to build" )
mark_as_advanced( ACTS_VECMEM_GIT_REPOSITORY ACTS_VECMEM_GIT_TAG )
FetchContent_Declare( VecMem
   GIT_REPOSITORY "${ACTS_VECMEM_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_VECMEM_GIT_TAG}" )

# Options used in the build of VecMem.
set( VECMEM_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the VecMem unit tests" )
set( VECMEM_BUILD_BENCHMARKING FALSE CACHE BOOL
   "Turn off the build of VecMem benchmarking" )   

# Now set up its build.
FetchContent_MakeAvailable( VecMem )
