# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray CUDA benchmarks")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Build benchmark library for multiple algebra plugins to create correct
# template instantiations
# Currently vc and smatrix is not supported on device
set(algebra_plugins "array")
if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebra_plugins "eigen")
endif()

# Set up a benchmark library for CUDA
foreach(algebra ${algebra_plugins})
    add_library(
        detray_benchmark_cuda_${algebra}
        STATIC
        "propagation_benchmark.hpp"
        "propagation_benchmark.cu"
    )

    add_library(
        detray::benchmark_cuda_${algebra}
        ALIAS detray_benchmark_cuda_${algebra}
    )

    target_link_libraries(
        detray_benchmark_cuda_${algebra}
        PUBLIC
            vecmem::cuda
            detray::benchmarks
            detray::test_common
            detray::core_${algebra}
    )
endforeach()
