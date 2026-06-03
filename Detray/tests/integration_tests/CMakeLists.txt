# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up all of the "host" tests.
if(DETRAY_BUILD_HOST)
    message(STATUS "Building detray host integration tests")

    add_subdirectory(cpu)
    add_subdirectory(io)
endif()

# Set up all of the "device" tests.
if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL OR DETRAY_BUILD_HIP)
    add_subdirectory(device)
endif()
