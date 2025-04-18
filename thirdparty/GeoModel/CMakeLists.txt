# This file is part of the Acts project.
#
# Copyright (C) 2024 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building GeoModel as part of the ACTS project" )

set( GEOMODEL_VERSION "${_acts_geomodel_version}")

set( GEOMODEL_SETUP_JSON OFF CACHE BOOL "Skip setting up json completely" )


# Declare where to get VecMem from.
set( ACTS_GEOMODEL_GIT_REPOSITORY "https://gitlab.cern.ch/GeoModelDev/GeoModel"
   CACHE STRING "Git repository to take GeoModel from" )
set( ACTS_GEOMODEL_GIT_TAG "${GEOMODEL_VERSION}" CACHE STRING "Version of GeoModel to build" )
mark_as_advanced( ACTS_GEOMODEL_GIT_REPOSITORY ACTS_GEOMODEL_GIT_TAG )
FetchContent_Declare( geomodel
   GIT_REPOSITORY "${ACTS_GEOMODEL_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_GEOMODEL_GIT_TAG}"
   PATCH_COMMAND git am ${CMAKE_CURRENT_SOURCE_DIR}/0001-Add-option-to-skip-setting-up-json-completely.patch
)

# Now set up its build.
FetchContent_MakeAvailable( geomodel )
