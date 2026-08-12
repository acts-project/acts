# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include(traccc-functions)
include(CheckCXXCompilerFlag)

# Turn on the correct setting for the __cplusplus macro with MSVC.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    traccc_add_flag( CMAKE_CXX_FLAGS "/Zc:__cplusplus" )
endif()

# Note that the C++ warning flags are *not* set here. They come from the shared
# ActsWarningFlags.cmake module and are applied once, per-directory, from
# Traccc/CMakeLists.txt -- so that they apply identically whether traccc is
# built standalone or as part of the Acts monorepo, and so that they never reach
# the dependencies set up in Traccc/extern.

# Set architecture level to v2 - enable generation of code using SSE4.2.
# But only use the flag if it's available from the compiler. And the user
# didn't specify some other -march flag yet.
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
    set(X86_FLAG "-march=x86-64-v2")
    check_cxx_compiler_flag("${X86_FLAG}" TRACCC_SSE42_SUPPORTED)
    if(TRACCC_SSE42_SUPPORTED AND (NOT "${CMAKE_CXX_FLAGS}" MATCHES "-march="))
        traccc_add_flag( CMAKE_CXX_FLAGS "${X86_FLAG}" )
    endif()
    unset(X86_FLAG)
endif()
