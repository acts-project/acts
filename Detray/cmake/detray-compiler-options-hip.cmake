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
        detray_add_flag( CMAKE_HIP_FLAGS_DEBUG "-G" )
        detray_add_flag( CMAKE_HIP_FLAGS "--expt-relaxed-constexpr" )
    endif()
endif()

# Warnings and warnings-as-errors come from the shared module. CMake picks the
# right flag for the platform in use (`-Werror` on amd, `-Werror all-warnings`
# on nvidia).
acts_apply_warning_flags(PROFILE detray LANGUAGES HIP)
