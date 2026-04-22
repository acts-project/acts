# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'detray::test_device' component")

# Set up a device test library
add_library(
    detray_test_device
    INTERFACE
    "device_fixture.hpp"
    "execute_host_test.hpp"
    "matrix_fixture.hpp"
    "propagator_test.hpp"
    "transform_fixture.hpp"
    "vector_fixture.hpp"
)
add_library(detray::test_device ALIAS detray_test_device)

target_link_libraries(
    detray_test_device
    INTERFACE
        detray::core_array
        detray::test_framework
        detray::test_utils
        detray::validation_utils
)

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()

if(DETRAY_BUILD_SYCL)
    add_subdirectory(sycl)
endif()
