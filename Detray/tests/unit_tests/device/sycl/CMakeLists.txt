# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray SYCL unit tests")

# Enable SYCL as a language.
vecmem_check_language( SYCL )

enable_language(SYCL)

# Set the default SYCL compiler flags.
include(detray-compiler-options-sycl)

# Define unit tests for multiple algebra plugins
# Currently, only 'array' is supported in SYCL code
set(algebras "array")

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_unit_test(sycl_${algebra}
       "algebra.sycl"
       LINK_LIBRARIES GTest::gtest_main vecmem::sycl detray::test_common detray::test_sycl
    )

    target_compile_definitions(
        detray_unit_test_sycl_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
