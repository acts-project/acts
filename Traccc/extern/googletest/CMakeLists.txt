# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required(VERSION 3.25)
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building GoogleTest as part of the TRACCC project")

# Disregard the global C++ flags, to avoid any unhelpful warnings.
set(CMAKE_CXX_FLAGS)

# Declare where to get GoogleTest from.
set(TRACCC_GOOGLETEST_SOURCE
    "URL;https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz;URL_MD5;7e11f6cfcf6498324ac82d567dcb891e"
    CACHE STRING
    "Source for GoogleTest, when built as part of this project"
)
mark_as_advanced(TRACCC_GOOGLETEST_SOURCE)
FetchContent_Declare(GoogleTest SYSTEM ${TRACCC_GOOGLETEST_SOURCE})

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
FetchContent_MakeAvailable(GoogleTest)

# Set up aliases for the GTest targets with the same name that they have
# when we find GTest pre-installed.
if(NOT TARGET GTest::gtest)
    add_library(GTest::gtest ALIAS gtest)
endif()
if(NOT TARGET GTest::gtest_main)
    add_library(GTest::gtest_main ALIAS gtest_main)
endif()
