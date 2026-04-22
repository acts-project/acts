# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray SYCL integration tests")

# Set the SYCL build flags.
include(detray-compiler-options-sycl)

# Enable SYCL as a language.
enable_language(SYCL)

# make unit tests for multiple algebras
# Currently eigen, vc and smatrix is not supported
set(algebras "array")

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_integration_test(sycl_${algebra}
       "propagator_sycl_kernel.hpp"
       "propagator.sycl"
       "propagator_kernel.sycl"
       LINK_LIBRARIES GTest::gtest_main vecmem::sycl detray::test_device
       detray::core detray::algebra_${algebra} detray::test_utils
    )
    target_compile_definitions(
        detray_integration_test_sycl_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
