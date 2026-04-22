# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building the Eigen3 algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_eigen algebra_eigen
   "include/algebra/eigen.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(detray_algebra_eigen INTERFACE Eigen3::Eigen detray::core)

target_compile_definitions(
    detray_algebra_eigen
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_EIGEN
)

detray_test_public_headers( detray_algebra_eigen
   "algebra/eigen.hpp"
)
