# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray CUDA unit tests")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Independent of algebra-plugin
detray_add_unit_test( cuda
   "container_cuda_kernel.hpp"
   "container_cuda.cpp"
   "container_cuda_kernel.cu"
   "utils_ranges_cuda_kernel.hpp"
   "utils_ranges_cuda.cpp"
   "utils_ranges_cuda_kernel.cu"
   LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray::core
)

# Define unit tests for multiple algebra plugins
# Currently, only 'array' and 'eigen' are supported in device code
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_unit_test(cuda_${algebra}
       "algebra.cu"
       "detector_cuda_kernel.hpp"
       "detector_cuda.cpp"
       "detector_cuda_kernel.cu"
       "mask_store_cuda_kernel.hpp"
       "mask_store_cuda.cpp"
       "mask_store_cuda_kernel.cu"
       "navigator_cuda_kernel.hpp"
       "navigator_cuda.cpp"
       "navigator_cuda_kernel.cu"
       "sf_finders_grid_cuda_kernel.hpp"
       "sf_finders_grid_cuda.cpp"
       "sf_finders_grid_cuda_kernel.cu"
       "transform_store_cuda_kernel.hpp"
       "transform_store_cuda.cpp"
       "transform_store_cuda_kernel.cu"
       LINK_LIBRARIES GTest::gtest_main
        CUDA::cudart vecmem::cuda covfie::cuda detray::core
       detray::algebra_${algebra} detray::test_common detray::test_utils
    )

    target_compile_definitions(
        detray_unit_test_cuda_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
