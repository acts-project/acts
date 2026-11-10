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

    # Generate debug symbols for the device code as well in a debug build.
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
    )
        detray_add_flag(CMAKE_HIP_FLAGS_DEBUG "-G")
        detray_add_flag(CMAKE_HIP_FLAGS "--expt-relaxed-constexpr")
    endif()
endif()

acts_apply_warning_flags(HIP ${DETRAY_FAIL_ON_WARNINGS})
