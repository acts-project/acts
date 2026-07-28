# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.21)

# Only set these compiler flags if we are the top level project.
if(PROJECT_IS_TOP_LEVEL)
    # Include the helper function(s).
    include(detray-functions)

    # Figure out the properties of CUDA being used.
    find_package(CUDAToolkit REQUIRED)

    # Turn on the correct setting for the __cplusplus macro with MSVC.
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
        detray_add_flag( CMAKE_CUDA_FLAGS "-Xcompiler /Zc:__cplusplus" )
    endif()

    # Set the CUDA architecture to build code for.
    set(CMAKE_CUDA_ARCHITECTURES
        "52"
        CACHE STRING
        "CUDA architectures to build device code for"
    )

    if("${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA")
        # Allow to use functions in device code that are constexpr, even if they are
        # not marked with __device__.
        detray_add_flag( CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr" )
    endif()

    # Make CUDA generate debug symbols for the device code as well in a debug
    # build.
    detray_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-G -src-in-ptx" )
    detray_add_flag( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-lineinfo -src-in-ptx" )
endif()

# Warnings and warnings-as-errors come from the shared module, and are applied
# in the directory that includes this file -- regardless of whether detray is
# the top level project. CMake itself picks the right flag for the compiler in
# use (`-Werror all-warnings` for nvcc >= 10.2, `-Werror` for clang-cuda).
acts_apply_warning_flags(PROFILE detray LANGUAGES CUDA)
