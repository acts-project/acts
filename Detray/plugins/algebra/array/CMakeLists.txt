# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building the std::array algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_array algebra_array
   "include/algebra/array.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(detray_algebra_array INTERFACE detray::core)

target_compile_definitions(
    detray_algebra_array
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_ARRAY
)

detray_test_public_headers( detray_algebra_array
   "algebra/array.hpp"
)
