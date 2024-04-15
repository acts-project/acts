# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building Thrust as part of the ACTS project" )

set( THRUST_VERSION "${_acts_thrust_version}")

# Declare where to get Thrust from.
set( ACTS_THRUST_GIT_REPOSITORY "https://github.com/NVIDIA/thrust.git"
   CACHE STRING "Git repository to take trust from" )
set( ACTS_THRUST_GIT_TAG "${THRUST_VERSION}" CACHE STRING "Version of thrust to build" )
mark_as_advanced( ACTS_THRUST_GIT_REPOSITORY ACTS_THRUST_GIT_TAG )
FetchContent_Declare( Thrust
   GIT_REPOSITORY "${ACTS_THRUST_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_THRUST_GIT_TAG}" )

# Options used in the build of Thrust.
set( THRUST_ENABLE_INSTALL_RULES TRUE CACHE BOOL
   "Install Thrust together with this project" )

# Get it into the current directory.
FetchContent_MakeAvailable( Thrust )

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Tell the user what's happening.
message( STATUS "Building Thrust as part of the Acts project" )