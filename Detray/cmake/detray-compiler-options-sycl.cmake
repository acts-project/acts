# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.21)

# Only set these compiler flags if we are the top level project.
if(PROJECT_IS_TOP_LEVEL)
    # Include the helper function(s).
    include(detray-functions)

    # Basic flags for all build modes.
    detray_add_flag( CMAKE_SYCL_FLAGS "-Wall" )
    detray_add_flag( CMAKE_SYCL_FLAGS "-Wextra" )
    detray_add_flag( CMAKE_SYCL_FLAGS "-Wno-unknown-cuda-version" )
    detray_add_flag( CMAKE_SYCL_FLAGS "-Wshadow" )
    detray_add_flag( CMAKE_SYCL_FLAGS "-Wunused-local-typedefs" )
    if(NOT WIN32)
        detray_add_flag( CMAKE_SYCL_FLAGS "-pedantic" )
    endif()

    # Avoid issues coming from MSVC<->DPC++ argument differences.
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
        detray_add_flag(CMAKE_SYCL_FLAGS
          "-Wno-unused-command-line-argument"
        )
    endif()

    # Fail on warnings, if asked for that behaviour.
    if(DETRAY_FAIL_ON_WARNINGS)
        detray_add_flag( CMAKE_SYCL_FLAGS "-Werror" )
    endif()
endif()
