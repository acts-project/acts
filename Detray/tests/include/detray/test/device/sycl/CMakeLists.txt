# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'detray::test_sycl' component")

# Enable SYCL as a language.
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
