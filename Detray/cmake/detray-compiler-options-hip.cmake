# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.21)

if(PROJECT_IS_TOP_LEVEL)
    include(detray-functions)

    #find HIP
    find_package(HIP REQUIRED)

    set(CMAKE_HIP_ARCHITECTURES gfx1031)

    # Basic flags for all build modes.
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    )
        detray_add_flag( CMAKE_HIP_FLAGS "-Wall" )
        detray_add_flag( CMAKE_HIP_FLAGS "-Wextra" )
        detray_add_flag( CMAKE_HIP_FLAGS "-Wshadow" )
        detray_add_flag( CMAKE_HIP_FLAGS "-Wunused-local-typedefs" )
        detray_add_flag( CMAKE_HIP_FLAGS "-pedantic" )
    endif()
    # Generate debug symbols for the device code as well in a debug build.
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
    )
        detray_add_flag( CMAKE_HIP_FLAGS_DEBUG "-G" )
        detray_add_flag( CMAKE_HIP_FLAGS "--expt-relaxed-constexpr" )
    endif()

    # Fail on warnings, if asked for that behaviour.
    if(DETRAY_FAIL_ON_WARNINGS)
        if(
            ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
            OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
        )
            detray_add_flag( CMAKE_HIP_FLAGS "-Werror" )
        elseif(
            ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvcc")
            OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
        )
            detray_add_flag( CMAKE_HIP_FLAGS "-Werror all-warnings" )
        endif()
    endif()
endif()
