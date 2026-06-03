# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.11)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

if(TARGET nlohmann_json::nlohmann_json)
    message(
        FATAL_ERROR
        "Setup of nlohmann JSON requested, but target `nlohmann_json::nlohmann_json` exists!"
    )
endif()

# Tell the user what's happening.
message(STATUS "Building nlohmann_json as part of the detray project")

# Declare where to get Eigen3 from.
set(DETRAY_NLOHMANN_SOURCE
    "URL;https://github.com/nlohmann/json/archive/refs/tags/v${DETRAY_NLOHMANN_VERSION}.tar.gz;URL_HASH;SHA1=${DETRAY_NLOHMANN_DIGEST}"
    CACHE STRING
    "Source for Nlohmann JSON, when built as part of this project"
)
mark_as_advanced(DETRAY_EIGEN_SOURCE)

# Mark the import as a system library on modern CMake versions
# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_NLOHMANN_JSON_SOURCE_FULL "${DETRAY_NLOHMANN_SOURCE};SYSTEM")
else()
    set(DETRAY_NLOHMANN_JSON_SOURCE_FULL "${DETRAY_NLOHMANN_SOURCE}")
endif()
mark_as_advanced(DETRAY_NLOHMANN_JSON_SOURCE_FULL)

FetchContent_Declare(nlohmann_json ${DETRAY_NLOHMANN_JSON_SOURCE_FULL})

# Now set up its build.
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install ON CACHE INTERNAL "")
FetchContent_MakeAvailable(nlohmann_json)
