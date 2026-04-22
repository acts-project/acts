# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::test_cpu' component")

# Set the CPU build flags.
include(detray-compiler-options-cpp)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    detray_test_cpu
    INTERFACE
    "algebra_fixture.hpp"
    "detector_consistency.hpp"
    "detector_scan.hpp"
    "material_scan.hpp"
    "navigation_validation.hpp"
    "navigation_validation.hpp"
    "toy_detector_test.hpp"
)

add_library(detray::test_cpu ALIAS detray_test_cpu)

target_link_libraries(
    detray_test_cpu
    INTERFACE
        vecmem::core
        detray::core_array
        detray::test_common
        detray::validation_utils
)
