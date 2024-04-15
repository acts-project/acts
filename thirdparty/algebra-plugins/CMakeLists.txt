# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Tell the user what's happening.
message( STATUS "Building algebra-plugins as part of the Acts project" )

set( ALGEBRAPLUGINS_VERSION "v${_acts_algebraplugins_version}")

# Declare where to get AlgebraPlugins from.
set( ACTS_ALGEBRAPLUGINS_GIT_REPOSITORY "https://github.com/acts-project/algebra-plugins.git"
   CACHE STRING "Git repository to take algebra-plugins from" )
set( ACTS_ALGEBRAPLUGINS_GIT_TAG "${ALGEBRAPLUGINS_VERSION}" CACHE STRING "Version of algebra-plugins to build" )
mark_as_advanced( ACTS_ALGEBRAPLUGINS_GIT_REPOSITORY ACTS_ALGEBRAPLUGINS_GIT_TAG )
FetchContent_Declare( AlgebraPlugins
   GIT_REPOSITORY "${ACTS_ALGEBRAPLUGINS_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_ALGEBRAPLUGINS_GIT_TAG}" )



# Options used in the build of Algebra Plugins.
set(ALGEBRA_PLUGINS_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the Algebra Plugins unit tests")
set(ALGEBRA_PLUGINS_INCLUDE_EIGEN TRUE CACHE BOOL
   "Turn on the build of algebra::eigen")
set(ALGEBRA_PLUGINS_INCLUDE_VC TRUE CACHE BOOL
   "Turn on the build of algebra::vc_array")
set(ALGEBRA_PLUGINS_INCLUDE_VECMEM TRUE CACHE BOOL
   "Turn on the build of algebra::vecmem_array")

# Build smatrix plugin only when it needs to be
# ROOT is not provided in native github compilers :/
set(ALGEBRA_PLUGINS_INCLUDE_SMATRIX ${DETRAY_SMATRIX_PLUGIN}
   CACHE BOOL "Turn on/off the build of algebra::smatrix")

set(ALGEBRA_PLUGINS_SETUP_EIGEN3 OFF CACHE BOOL
   "Do not have Algebra Plugins set up Eigen3 for itself")
set(ALGEBRA_PLUGINS_SETUP_VC ON CACHE BOOL
   "Have Algebra Plugins set up Vc for itself")
set(ALGEBRA_PLUGINS_USE_SYSTEM_VC OFF CACHE BOOL
   "Have Algebra Plugins build Vc itself")
set(ALGEBRA_PLUGINS_SETUP_VECMEM FALSE CACHE BOOL
   "Do not set up VecMem in Algebra Plugins")
set(ALGEBRA_PLUGINS_SETUP_GOOGLETEST FALSE CACHE BOOL
   "Do not set up GoogleTest in Algebra Plugins")

set(ALGEBRA_PLUGINS_SETUP_BENCHMARK FALSE CACHE BOOL "Do not setup Algebra Plugins benchmark")

# Now set up its build.
FetchContent_MakeAvailable( AlgebraPlugins )