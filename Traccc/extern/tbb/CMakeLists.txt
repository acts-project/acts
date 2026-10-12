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
message(STATUS "Building TBB as part of the TRACCC project")

# Disregard the global C++ flags, to avoid any unhelpful warnings.
set(CMAKE_CXX_FLAGS)

# Declare where to get TBB from.
set(TRACCC_TBB_SOURCE
    "URL;https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2022.0.0.tar.gz;URL_MD5;78ec44cecf3cd78c4984e61f1bb93134"
    CACHE STRING
    "Source for TBB, when built as part of this project"
)
mark_as_advanced(TRACCC_TBB_SOURCE)
FetchContent_Declare(TBB SYSTEM ${TRACCC_TBB_SOURCE} OVERRIDE_FIND_PACKAGE)

# Options used in the build of TBB.
set(TBB_TEST FALSE CACHE BOOL "Turn off the TBB tests")
set(TBB_STRICT FALSE CACHE BOOL "Do not throw errors on compiler warnings")

# Make TBB work without warnings with modern CMake versions.
set(CMAKE_POLICY_VERSION_MINIMUM "3.10")

# Get it into the current directory.
FetchContent_MakeAvailable(TBB)
