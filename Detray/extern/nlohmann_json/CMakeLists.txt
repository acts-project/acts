# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.11)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
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
