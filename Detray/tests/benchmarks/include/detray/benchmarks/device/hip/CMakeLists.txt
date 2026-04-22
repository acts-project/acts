# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.21) # HIP language support requires minimum 3.21

find_package(HIPToolkit)

# Enable HIP as a language.
enable_language(HIP)

# Set the HIP build flags.
include(detray-compiler-options-hip)

# Build benchmark library for multiple algebra plugins to create correct
# template instantiations
# Currently vc and smatrix is not supported on device
set(algebra_plugins "array")
if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebra_plugins "eigen")
endif()

# Set up a benchmark library for HIP
foreach(algebra ${algebra_plugins})
    add_library(
        detray_benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}
        STATIC
        "propagation_benchmark.hpp"
        "propagation_benchmark.hip"
    )

    add_library(
        detray::benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}
        ALIAS detray_benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}
    )

    target_link_libraries(
        detray_benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}
        PUBLIC
            HIP::hiprt
            vecmem::hip
            detray::benchmarks
            detray::test_common
            detray::core_${algebra}
    )
endforeach()
