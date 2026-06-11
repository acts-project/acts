# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::test_sycl' component")

# Enable SYCL as a language.
include(vecmem-check-language) # Needed for the SYCL language.
enable_language(SYCL)

# Set the default SYCL compiler flags.
include(detray-compiler-options-sycl)

# Helper library for the SYCL tests.
add_library(
    detray_test_sycl
    INTERFACE
    "algebra_test_suite.hpp"
    "execute_sycl_test.hpp"
    "sycl_test_fixture.hpp"
)

add_library(detray::test_sycl ALIAS detray_test_sycl)

target_link_libraries(
    detray_test_sycl
    INTERFACE GTest::gtest vecmem::sycl detray::core_array detray::test_device
)
