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
message(STATUS "Building Alpaka as part of the TRACCC project")

# Disregard the global C++ flags, to avoid any unhelpful warnings.
set(CMAKE_CXX_FLAGS)

# Declare where to get Alpaka from.
set(TRACCC_ALPAKA_SOURCE
    "URL;https://github.com/alpaka-group/alpaka/archive/refs/tags/2.2.0.tar.gz;URL_MD5;c12ee1495974ec911c1a43e0ff4d8925"
    CACHE STRING
    "Source for Alpaka, when built as part of this project"
)

set(TRACCC_ALPAKA_PATCH
    ""
    CACHE STRING
    "Patch(es) for alpaka, when built as part of this project"
)

mark_as_advanced(TRACCC_ALPAKA_SOURCE TRACCC_ALPAKA_PATCH)
FetchContent_Declare(
    Alpaka
    SYSTEM
    ${TRACCC_ALPAKA_SOURCE}
    ${TRACCC_ALPAKA_PATCH}
)

# Turn off build testing.
set(BUILD_TESTING FALSE)

# Get it into the current directory.
FetchContent_MakeAvailable(Alpaka)
