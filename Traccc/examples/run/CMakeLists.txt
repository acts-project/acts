# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Project include(s).
include(traccc-compiler-options-cpp)

# # Add all the subdirectories that can be built.
add_subdirectory(common)
add_subdirectory(validation)

if(TRACCC_BUILD_CUDA)
    add_subdirectory(cuda)
endif()

if(TRACCC_BUILD_HIP)
    add_subdirectory(hip)
endif()

if(TRACCC_BUILD_SYCL)
    add_subdirectory(sycl)
endif()

if(TRACCC_BUILD_ALPAKA)
    add_subdirectory(alpaka)
endif()
