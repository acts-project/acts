# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
