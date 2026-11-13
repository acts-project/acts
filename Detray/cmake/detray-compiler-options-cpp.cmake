# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.21)

# Only set these compiler flags if we are the top level project.
if(PROJECT_IS_TOP_LEVEL)
    # Include the helper function(s).
    include(detray-functions)

    # Turn on the correct setting for the __cplusplus macro with MSVC.
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
        detray_add_flag(CMAKE_CXX_FLAGS "/Zc:__cplusplus")
    endif()

    # Respect infinity expressions for IntelLLVM
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "IntelLLVM")
        detray_add_flag(CMAKE_CXX_FLAGS "-fhonor-infinities")
    endif()
endif()

# C++ warnings are applied in Detray/CMakeLists.txt after dependency setup.
