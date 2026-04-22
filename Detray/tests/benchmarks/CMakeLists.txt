# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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
