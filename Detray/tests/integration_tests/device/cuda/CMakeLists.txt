# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray CUDA integration tests")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_integration_test(cuda_${algebra}
       "propagator_cuda_kernel.hpp"
       "propagator_cuda.cpp"
       "propagator_cuda_kernel.cu"
       LINK_LIBRARIES GTest::gtest_main CUDA::cudart vecmem::cuda covfie::cuda
                       detray::core detray::algebra_${algebra} detray::test_device
                       detray::test_common
    )

    target_compile_definitions(
        detray_integration_test_cuda_${algebra}
        PRIVATE ${algebra}=${algebra}
    )

    set_tests_properties(
        detray_integration_test_cuda_${algebra}
        PROPERTIES
            DEPENDS "detray_unit_test_cuda;detray_unit_test_cuda_${algebra}"
    )
endforeach()

# CUDA navigation validation for the telescope detector
detray_add_integration_test(telescope_detector_cuda
                           "telescope_navigation_validation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            vecmem::cuda detray::core_array
                            detray::test_cuda detray::test_common
)

set_tests_properties(
    detray_integration_test_telescope_detector_cuda
    PROPERTIES
        DEPENDS
            "detray_integration_test_cuda_array,detray_integration_test_telescope_detector"
)

# CUDA navigation validation for the toy detector
detray_add_integration_test(toy_detector_cuda
                           "toy_detector_navigation_validation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            vecmem::cuda covfie::cuda detray::core_array
                            detray::test_cuda detray::test_common
)

set_tests_properties(
    detray_integration_test_toy_detector_cuda
    PROPERTIES
        DEPENDS
            "detray_integration_test_cuda_array,detray_integration_test_toy_detector"
)

# CUDA navigation validation for the wire chamber
detray_add_integration_test(wire_chamber_cuda
                           "wire_chamber_navigation_validation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            vecmem::cuda detray::core_array
                            detray::test_cuda detray::test_common
)

set_tests_properties(
    detray_integration_test_wire_chamber_cuda
    PROPERTIES
        DEPENDS
            "detray_integration_test_cuda_array,detray_integration_test_wire_chamber"
)
