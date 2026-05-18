# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building the SMatrix algebra plugin")

find_package(ROOT COMPONENTS Core MathCore Smatrix REQUIRED)

# Set up the library.
detray_add_library( detray_algebra_smatrix algebra_smatrix
   "include/algebra/smatrix.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
)

target_link_libraries(
    detray_algebra_smatrix
    INTERFACE ROOT::Core ROOT::MathCore ROOT::Smatrix detray::core
)

target_compile_definitions(
    detray_algebra_smatrix
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_SMATRIX
)

detray_test_public_headers( detray_algebra_smatrix
   "algebra/smatrix.hpp"
)
