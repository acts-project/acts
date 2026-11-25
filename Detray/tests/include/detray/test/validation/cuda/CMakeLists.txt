# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::validation_cuda' component")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Set up a validation library for CUDA
add_library(
    detray_validation_cuda
    STATIC
    "material_recorder.hpp"
    "material_recorder.cu"
)

add_library(detray::validation_cuda ALIAS detray_validation_cuda)

target_link_libraries(
    detray_validation_cuda
    PUBLIC
        CUDA::cudart
        vecmem::cuda
        covfie::cuda
        detray::core_array
        detray::test_device
        detray::validation
)
