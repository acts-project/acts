# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set the common C++ flags.
include(detray-compiler-options-cpp)

include_directories(
    SYSTEM
    $<TARGET_PROPERTY:covfie::core,INTERFACE_INCLUDE_DIRECTORIES>
)

# Set up a common benchmark library.
file(
    GLOB _detray_benchmarks_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/benchmarks/*.hpp"
)

add_library(detray_benchmarks INTERFACE "${_detray_benchmarks_headers}")
add_library(detray::benchmarks ALIAS detray_benchmarks)

target_include_directories(
    detray_benchmarks
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    detray_benchmarks
    INTERFACE benchmark::benchmark vecmem::core detray::core detray::test_common
)

unset(_detray_benchmarks_headers)

# Set up the host/cpu benchmarks.
if(DETRAY_BUILD_HOST)
    add_subdirectory(cpu)
    add_subdirectory(include/detray/benchmarks/cpu)
endif()

# Set up all of the "device" benchmarks.

if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_HIP)
    add_subdirectory(include/detray/benchmarks/device)
endif()

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()
