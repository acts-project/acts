# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Include the helper function(s).
include(traccc-functions)

# Only tweak the flags for the Intel compiler.
if(
    NOT (
        ("${CMAKE_SYCL_COMPILER_ID}" STREQUAL "IntelLLVM")
        OR ("${CMAKE_SYCL_COMPILER_ID}" MATCHES "Clang")
    )
)
    return()
endif()

acts_apply_warning_flags(SYCL ${TRACCC_FAIL_ON_WARNINGS})
