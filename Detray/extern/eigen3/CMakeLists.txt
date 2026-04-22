# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Eigen3 as part of the Detray project")

# Declare where to get Eigen3 from.
set(DETRAY_EIGEN_SOURCE
    "URL;https://gitlab.com/libeigen/eigen/-/archive/${DETRAY_EIGEN_VERSION}/eigen-${DETRAY_EIGEN_VERSION}.tar.bz2;URL_MD5;${DETRAY_EIGEN_DIGEST}"
    CACHE STRING
    "Source for Eigen, when built as part of this project"
)
mark_as_advanced(DETRAY_EIGEN_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_EIGEN_SOURCE_FULL "${DETRAY_EIGEN_SOURCE};SYSTEM")
else()
    set(DETRAY_EIGEN_SOURCE_FULL "${DETRAY_EIGEN_SOURCE}")
endif()
mark_as_advanced(DETRAY_EIGEN_SOURCE_FULL)

FetchContent_Declare(Eigen3 ${DETRAY_EIGEN_SOURCE_FULL})

# Turn off the unit tests for Eigen3.
if(DEFINED CACHE{BUILD_TESTING})
    set(_buildTestingValue ${BUILD_TESTING})
endif()
set(BUILD_TESTING FALSE CACHE INTERNAL "Forceful setting of BUILD_TESTING")

# Get it into the current directory.
FetchContent_MakeAvailable(Eigen3)

# Reset the BUILD_TESTING variable.
if(DEFINED _buildTestingValue)
    set(BUILD_TESTING
        ${_buildTestingValue}
        CACHE BOOL
        "Turn tests on/off"
        FORCE
    )
    unset(_buildTestingValue)
else()
    unset(BUILD_TESTING CACHE)
endif()

# Treat the Eigen headers as "system headers", to avoid getting warnings from
# them.
get_target_property(_incDirs eigen INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(eigen SYSTEM INTERFACE ${_incDirs})
