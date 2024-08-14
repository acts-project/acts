# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message( STATUS "Building tracc as part of the Acts project" )

set( TRACCC_VERSION "${_acts_traccc_version}")

# Declare where to get traccc from.
file(
    DOWNLOAD
        https://github.com/acts-project/traccc/commit/d2a8754242ddfa84f03fcab837205f548279e4ad.diff
        ${CMAKE_CURRENT_BINARY_DIR}/traccc-patch.diff
    EXPECTED_HASH
        SHA256=ae9140fbfc51b28bbd5857b59bec8af36d09d50d4fd3756627b5ab76c2014068
)

# Declare where to get traccc from.
FetchContent_Declare(
    traccc
    ${ACTS_TRACCC_SOURCE}
    PATCH_COMMAND patch -p1 -i ${CMAKE_CURRENT_BINARY_DIR}/traccc-patch.diff
)

if(ACTS_CUSTOM_SCALARTYPE)
   set(ACTS_TRACCC_SCALARTYPE ${ACTS_CUSTOM_SCALARTYPE})
else()
   set(ACTS_TRACCC_SCALARTYPE "double")
endif()

set( TRACCC_CUSTOM_SCALARTYPE "${ACTS_TRACCC_SCALARTYPE}" CACHE STRING
   "Scalar type to use in the traccc code" )

set( TRACCC_SETUP_VECMEM OFF CACHE BOOL
   "Do not set up Actsvg as part of Traccc" )
set( TRACCC_SETUP_EIGEN3 OFF CACHE BOOL
   "Do not set up Eigen3 as part of Traccc" )
set( TRACCC_SETUP_THRUST ON CACHE BOOL
    "Do not set up Thrust as part of Traccc" )
set( TRACCC_SETUP_ALGEBRA_PLUGINS OFF CACHE BOOL
    "Do not set up Algebra Plugins as part of Traccc" )
set( TRACCC_SETUP_COVFIE OFF CACHE BOOL
    "Do not set up Covfie as part of Traccc" )
set( TRACCC_SETUP_DFELIBS OFF CACHE BOOL
    "Do not set up dfelibs as part of Traccc" )
set( TRACCC_SETUP_DETRAY ON CACHE BOOL
    "Set up Detray as part of Traccc" )
set( TRACCC_SETUP_ACTS OFF CACHE BOOL
    "Do not set up ACTS as part of Traccc" )
set( TRACCC_SETUP_TBB OFF CACHE BOOL
    "Do not set up TBB as part of Traccc" )

set( TRACCC_BUILD_TESTING OFF CACHE BOOL
    "Turn off the build of the Traccc unit tests" )
set( TRACCC_BUILD_EXAMPLES OFF CACHE BOOL
    "Turn off the build of the Traccc examples" )

# Now set up its build.
FetchContent_MakeAvailable( traccc )
