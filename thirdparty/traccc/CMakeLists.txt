# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building traccc as part of the Acts project")

set(TRACCC_VERSION "${_acts_traccc_version}")

# Declare where to get traccc from.
FetchContent_Declare(
    traccc
    SYSTEM
    ${ACTS_TRACCC_SOURCE}
    PATCH_COMMAND
        patch -p1 <
        ${CMAKE_CURRENT_SOURCE_DIR}/0001-Fix-kalman-fitter-state-ref-order.patch
)

set(TRACCC_CUSTOM_SCALARTYPE
    "float"
    CACHE STRING
    "Scalar type to use in the traccc code"
)

set(TRACCC_SETUP_VECMEM OFF CACHE BOOL "Do not set up VecMem as part of Traccc")
set(TRACCC_SETUP_EIGEN3 OFF CACHE BOOL "Do not set up Eigen3 as part of Traccc")
set(TRACCC_SETUP_THRUST ON CACHE BOOL "Do not set up Thrust as part of Traccc")

set(TRACCC_SETUP_COVFIE OFF CACHE BOOL "Do not set up Covfie as part of Traccc")
set(TRACCC_SETUP_DFELIBS
    ON
    CACHE BOOL
    "Do not set up dfelibs as part of Traccc"
)
set(TRACCC_SETUP_DETRAY OFF CACHE BOOL "Set up Detray as part of Traccc")
set(TRACCC_SETUP_ACTS OFF CACHE BOOL "Do not set up ACTS as part of Traccc")

set(TRACCC_BUILD_CUDA
    ${ACTS_ENABLE_CUDA}
    CACHE BOOL
    "Build Traccc with CUDA support"
)

set(TRACCC_BUILD_IO
    ON
    CACHE BOOL
    "Turn off the build of the Traccc IO code THIS IS TEMPORARY"
)

set(TRACCC_BUILD_PERFORMANCE
    OFF
    CACHE BOOL
    "Turn off the build of the Traccc performance code"
)

set(TRACCC_BUILD_SIMULATION
    OFF
    CACHE BOOL
    "Turn off the build of the Traccc simulation code"
)

set(TRACCC_BUILD_TESTING
    OFF
    CACHE BOOL
    "Turn off the build of the Traccc unit tests"
)

set(TRACCC_SETUP_TBB OFF CACHE BOOL "Do not set up TBB as part of Traccc")

set(TRACCC_BUILD_EXAMPLES
    OFF
    CACHE BOOL
    "Turn off the build of the Traccc examples"
)

set(TRACCC_BUILD_BENCHMARKS
    OFF
    CACHE BOOL
    "Turn off the build of the Traccc benchmarks"
)

# Now set up its build.
FetchContent_MakeAvailable(traccc)
