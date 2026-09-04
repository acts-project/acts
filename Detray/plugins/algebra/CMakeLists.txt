# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Temporary setting for the scalar type, until it can be removed.
set(DETRAY_CUSTOM_SCALARTYPE
    "float"
    CACHE STRING
    "Scalar type to use in the Detray code"
)

message(STATUS "Using detray scalar type: ${DETRAY_CUSTOM_SCALARTYPE}")

# Set up all enabled libraries.
if(DETRAY_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()

if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()

if(DETRAY_FASTOR_PLUGIN)
    add_subdirectory(fastor)
endif()

if(DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()

if(DETRAY_VC_AOS_PLUGIN)
    add_subdirectory(vc_aos)
endif()

if(DETRAY_VC_SOA_PLUGIN)
    add_subdirectory(vc_soa)
endif()
