# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.11)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Tell the user what's happening.
message(STATUS "Building Google Benchmark as part of the Detray project")

# Declare where to get Google Benchmark from.
set(DETRAY_BENCHMARK_SOURCE
    "URL;https://github.com/google/benchmark/archive/refs/tags/v${DETRAY_BENCHMARK_VERSION}.tar.gz;URL_MD5;${DETRAY_BENCHMARK_DIGEST}"
    CACHE STRING
    "Source for Google Benchmark, when built as part of this project"
)
mark_as_advanced(DETRAY_BENCHMARK_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_BENCHMARK_SOURCE_FULL "${DETRAY_BENCHMARK_SOURCE};SYSTEM")
else()
    set(DETRAY_BENCHMARK_SOURCE_FULL "${DETRAY_BENCHMARK_SOURCE}")
endif()
mark_as_advanced(DETRAY_BENCHMARK_SOURCE_FULL)

FetchContent_Declare(Benchmark ${DETRAY_BENCHMARK_SOURCE_FULL})

# Options used in the build of Google Benchmark.
set(BENCHMARK_ENABLE_TESTING
    OFF
    CACHE BOOL
    "Turn off the tests in Google Benchmark"
)
set(BENCHMARK_ENABLE_WERROR
    OFF
    CACHE BOOL
    "Turn off the -Werror for Release build"
)

# Get it into the current directory.
FetchContent_MakeAvailable(Benchmark)

# Set up an alias for the Google Benchmark target with the same name that it
# has when we find it pre-installed.
add_library(benchmark::benchmark ALIAS benchmark)
