# This file is part of the Acts project.
#
# Copyright (C) 2022 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building actsvg as part of the ACTS project" )

set( ACTS_VERSION "v${_acts_actsvg_version}")

# Declare where to get VecMem from.
set( ACTS_ACTSVG_GIT_REPOSITORY "https://github.com/acts-project/actsvg.git"
   CACHE STRING "Git repository to take actsvg from" )
set( ACTS_ACTSVG_GIT_TAG "${ACTS_VERSION}" CACHE STRING "Version of actsvg to build" )
mark_as_advanced( ACTS_ACTSVG_GIT_REPOSITORY ACTS_ACTSVG_GIT_TAG )
FetchContent_Declare( actsvg
   GIT_REPOSITORY "${ACTS_ACTSVG_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_ACTSVG_GIT_TAG}" )

# Now set up its build.
FetchContent_MakeAvailable( actsvg )
