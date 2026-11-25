# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set the common C++ flags.
include(detray-compiler-options-cpp)

include_directories(
    SYSTEM
    $<TARGET_PROPERTY:covfie::core,INTERFACE_INCLUDE_DIRECTORIES>
)

# Build extended validation tools (actsvg and csv dependencies)
message(STATUS "Building 'detray::validation' component")

if(NOT TARGET detray::svgtools)
    message(
        STATUS
        "detray::validation is missing dependencies: 'detray::svgtools' not available"
    )
    return()
endif()

# Check that ACTS_FORCE_ASSERTIONS works
acts_add_assertion_tests()

file(
    GLOB _detray_validation_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/test/validation/*.hpp"
    "include/detray/test/validation/cpu/*.hpp"
)

detray_add_library(detray_validation validation ${_detray_validation_headers})

target_link_libraries(
    detray_validation
    INTERFACE
        vecmem::core
        detray::core
        detray::test_utils
        detray::test_tools
        detray::csv_io
        detray::svgtools
)

unset(_detray_validation_headers)

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()
