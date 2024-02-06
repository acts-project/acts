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
message( STATUS "Building Detray as part of the Acts project" )

# Declare where to get Detray from.
set( ACTS_DETRAY_GIT_REPOSITORY "https://github.com/acts-project/detray.git"
   CACHE STRING "Git repository to take Detray from" )
set( ACTS_DETRAY_GIT_TAG "v0.51.0" CACHE STRING "Version of Detray to build" )
mark_as_advanced( ACTS_DETRAY_GIT_REPOSITORY ACTS_DETRAY_GIT_TAG )
FetchContent_Declare( Detray
   GIT_REPOSITORY "${ACTS_DETRAY_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_DETRAY_GIT_TAG}" )

# Now set up its build.
FetchContent_MakeAvailable(Detray)
