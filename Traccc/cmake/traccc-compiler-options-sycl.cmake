# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Include the helper function(s).
include(traccc-functions)

# Only tweak the flags for the Intel compiler.
if(
    NOT (
        ("${CMAKE_SYCL_COMPILER_ID}" STREQUAL "IntelLLVM")
        OR ("${CMAKE_SYCL_COMPILER_ID}" MATCHES "Clang")
    )
)
    return()
endif()

# Basic flags for all build modes.
foreach(mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wall")
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wextra")
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wno-unknown-cuda-version")
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wshadow")
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wunused-local-typedefs")
    traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Wconversion")
endforeach()

if(NOT WIN32)
    foreach(mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
        traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-pedantic")
    endforeach()
endif()

# Fail on warnings, if asked for that behaviour.
if(TRACCC_FAIL_ON_WARNINGS)
    foreach(mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
        traccc_add_flag(CMAKE_SYCL_FLAGS_${mode} "-Werror")
    endforeach()
endif()

# Avoid issues coming from MSVC<->DPC++ argument differences.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    foreach(mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
        traccc_add_flag(CMAKE_SYCL_FLAGS_${mode}
         "-Wno-unused-command-line-argument"
        )
    endforeach()
endif()
