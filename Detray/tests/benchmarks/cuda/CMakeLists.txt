# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Build benchmarks for multiple algebra plugins
# Currently vc and smatrix is not supported on device
set(algebra_plugins "array")
if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebra_plugins "eigen")
endif()

foreach(algebra ${algebra_plugins})
    message(STATUS "Building 'detray::benchmark_cuda_${algebra}' component")

    detray_add_executable(benchmark_cuda_propagation_${algebra}
      "propagation.cpp"
       LINK_LIBRARIES detray::benchmark_cuda_${algebra} detray::core_${algebra} vecmem::cuda detray::test_common
    )

    target_compile_definitions(
        detray_benchmark_cuda_propagation_${algebra}
        PRIVATE ${algebra}=${algebra}
    )

    target_compile_options(
        detray_benchmark_cuda_propagation_${algebra}
        PRIVATE "-march=native" "-ftree-vectorize"
    )
endforeach()
