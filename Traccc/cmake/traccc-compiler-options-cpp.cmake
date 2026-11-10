# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Include the helper function(s).
include(traccc-functions)
include(CheckCXXCompilerFlag)

# Turn on the correct setting for the __cplusplus macro with MSVC.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    traccc_add_flag(CMAKE_CXX_FLAGS "/Zc:__cplusplus")
endif()

acts_apply_warning_flags(CXX ${TRACCC_FAIL_ON_WARNINGS})

# Set architecture level to v2 - enable generation of code using SSE4.2.
# But only use the flag if it's available from the compiler. And the user
# didn't specify some other -march flag yet.
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
    set(X86_FLAG "-march=x86-64-v2")
    check_cxx_compiler_flag("${X86_FLAG}" TRACCC_SSE42_SUPPORTED)
    if(TRACCC_SSE42_SUPPORTED AND (NOT "${CMAKE_CXX_FLAGS}" MATCHES "-march="))
        traccc_add_flag(CMAKE_CXX_FLAGS "${X86_FLAG}")
    endif()
    unset(X86_FLAG)
endif()
