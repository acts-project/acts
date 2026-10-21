# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Tell the user what's happening.
message(STATUS "Building pybind11 as part of the Detray project")

# Declare where to get pybind11 from.
set(DETRAY_PYBIND11_SOURCE
    "GIT_REPOSITORY;https://github.com/pybind/pybind11.git;GIT_TAG;v${DETRAY_PYBIND11_VERSION}"
    CACHE STRING
    "Source for pybind11, when built as part of this project"
)
mark_as_advanced(DETRAY_PYBIND11_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_PYBIND11_SOURCE_FULL "${DETRAY_PYBIND11_SOURCE};SYSTEM")
else()
    set(DETRAY_PYBIND11_SOURCE_FULL "${DETRAY_PYBIND11_SOURCE}")
endif()
mark_as_advanced(DETRAY_PYBIND11_SOURCE_FULL)

# Declare where to get pybind from.
FetchContent_Declare(pybind11 ${DETRAY_PYBIND11_SOURCE_FULL})

# Set up its build.
set(PYBIND11_TEST OFF CACHE BOOL "Turn off the pybind11 test suite build")
set(PYBIND11_INSTALL OFF CACHE BOOL "Turn off the pybind11 installation")

# Get it into the current directory.
FetchContent_MakeAvailable(pybind11)
