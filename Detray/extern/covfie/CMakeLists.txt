# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Tell the user what's happening.
message(STATUS "Fetching covfie as part of the Detray project")

# Declare where to get covfie from.
set(DETRAY_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v${DETRAY_COVFIE_VERSION}.tar.gz;URL_MD5;${DETRAY_COVFIE_DIGEST}"
    CACHE STRING
    "Source for covfie, when built as part of this project"
)
mark_as_advanced(DETRAY_COVFIE_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_COVFIE_SOURCE_FULL "${DETRAY_COVFIE_SOURCE};SYSTEM")
else()
    set(DETRAY_COVFIE_SOURCE_FULL "${DETRAY_COVFIE_SOURCE}")
endif()
mark_as_advanced(DETRAY_COVFIE_SOURCE_FULL)

FetchContent_Declare(covfie ${DETRAY_COVFIE_SOURCE_FULL})

# Options used for covfie.
set(COVFIE_BUILD_EXAMPLES OFF CACHE BOOL "Build covfie examples")
set(COVFIE_BUILD_TESTS OFF CACHE BOOL "Build covfie tests")
set(COVFIE_BUILD_BENCHMARKS OFF CACHE BOOL "Build covfie benchmarks")

set(COVFIE_PLATFORM_CPU ON CACHE BOOL "Enable covfie CPU platform")

set(COVFIE_PLATFORM_HIP
    ${DETRAY_BUILD_HIP}
    CACHE BOOL
    "Enable covfie HIP platform"
)
set(COVFIE_PLATFORM_CUDA
    ${DETRAY_BUILD_CUDA}
    CACHE BOOL
    "Enable covfie CUDA platform"
)

set(COVFIE_REQUIRE_CXX20 OFF CACHE BOOL "Enable covfie C++20 requirement")
set(COVFIE_QUIET ON CACHE BOOL "Quiet covfie feature warnings")

# Get it into the current directory.
FetchContent_MakeAvailable(covfie)
