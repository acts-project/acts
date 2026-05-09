# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::test_cuda' component")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    detray_test_cuda
    STATIC
    "algebra_test_suite.cuh"
    "bfield.hpp"
    "cuda_test_fixture.hpp"
    "execute_cuda_test.cuh"
    "material_validation.hpp"
    "material_validation.cu"
    "navigation_validation.hpp"
    "navigation_validation.cu"
)

add_library(detray::test_cuda ALIAS detray_test_cuda)

target_link_libraries(
    detray_test_cuda
    PUBLIC
        CUDA::cudart
        vecmem::cuda
        covfie::cuda
        detray::core_array
        detray::test_device
        detray::test_cpu
        detray::validation_utils
)
