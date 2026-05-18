# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building detray CUDA tutorials")

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Examples using the array algebra plugin.
detray_add_tutorial( detector_construction_cuda
   "detector_construction.cpp"
   "detector_construction.hpp"
   "detector_construction.cu"
   LINK_LIBRARIES vecmem::cuda detray::core_array detray::tutorial
                  detray::test_common
)
target_compile_definitions(
    detray_tutorial_detector_construction_cuda
    PRIVATE array=array
)

detray_add_tutorial( propagation_cuda
   "propagation.cpp"
   "propagation.hpp"
   "propagation_kernel.cu"
   LINK_LIBRARIES CUDA::cudart vecmem::cuda detray::core_array covfie::core covfie::cuda
                  detray::tutorial detray::test_common detray::test_cuda
)
target_compile_definitions(detray_tutorial_propagation_cuda PRIVATE array=array)
