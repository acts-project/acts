# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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
