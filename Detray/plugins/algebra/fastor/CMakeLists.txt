# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building the Fastor algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_fastor algebra_fastor
   "include/algebra/fastor.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(
    detray_algebra_fastor
    INTERFACE Fastor::Fastor detray::core
)

target_compile_definitions(
    detray_algebra_fastor
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_FASTOR
)

# HACK: Fastor really requires FMA to be enabled, so we turn it on
# unconditionally. In Debug mode, it also requires at least -O1 and
# -fexpensive-optimizations.
target_compile_options(
    detray_algebra_fastor
    INTERFACE
        "-mfma"
        "$<$<CONFIG:Debug>:-O1>"
        "$<$<CONFIG:Debug>:-fexpensive-optimizations>"
)

detray_test_public_headers( detray_algebra_fastor
   "algebra/fastor.hpp"
)
