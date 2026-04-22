# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.11)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Tell the user what's happening.
message(STATUS "Building GoogleTest as part of the Detray project")

# Declare where to get GoogleTest from.
set(DETRAY_GOOGLETEST_SOURCE
    "URL;https://github.com/google/googletest/archive/refs/tags/v${DETRAY_GOOGLETEST_VERSION}.tar.gz;URL_MD5;${DETRAY_GOOGLETEST_DIGEST}"
    CACHE STRING
    "Source for GoogleTest, when built as part of this project"
)
mark_as_advanced(DETRAY_GOOGLETEST_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_GOOGLETEST_SOURCE_FULL "${DETRAY_GOOGLETEST_SOURCE};SYSTEM")
else()
    set(DETRAY_GOOGLETEST_SOURCE_FULL "${DETRAY_GOOGLETEST_SOURCE}")
endif()
mark_as_advanced(DETRAY_GOOGLETEST_SOURCE_FULL)

FetchContent_Declare(GoogleTest ${DETRAY_GOOGLETEST_SOURCE_FULL})

# Options used in the build of GoogleTest.
set(BUILD_GMOCK TRUE CACHE BOOL "Turn off the build of GMock")
set(INSTALL_GTEST FALSE CACHE BOOL "Turn off the installation of GoogleTest")
if(WIN32)
    set(gtest_force_shared_crt
        TRUE
        CACHE BOOL
        "Use shared (DLL) run-time library, even with static libraries"
    )
endif()

# Silence some warnings with modern versions of CMake on macOS.
set(CMAKE_MACOSX_RPATH TRUE)

# Get it into the current directory.
FetchContent_Populate(GoogleTest)
add_subdirectory(
    "${googletest_SOURCE_DIR}"
    "${googletest_BINARY_DIR}"
    EXCLUDE_FROM_ALL
)

# Set up aliases for the GTest targets with the same name that they have
# when we find GTest pre-installed.
add_library(GTest::gtest ALIAS gtest)
add_library(GTest::gtest_main ALIAS gtest_main)
