# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::detectors' component")

# Set the common C++ flags.
include(detray-compiler-options-cpp)

# Set up a detector library, which contains the different detector types that
# are currently available by default.
file(
    GLOB _detray_detector_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/detectors/*.hpp"
)

detray_add_library( detray_detectors detectors ${_detray_detectors_headers})
target_link_libraries(
    detray_detectors
    INTERFACE vecmem::core detray::core detray::io
)

unset(_detray_detector_headers)

# Compile all headers one-by-one.
detray_compile_headers(
    detectors
    detray::detectors
    "include/detray/detectors/*.hpp"
    "include/detray/detectors/**/*.hpp"
)
if(TARGET detray_detectors_HEADERS)
    # The detector headers need an algebra plugin to compile.
    target_link_libraries(
        detray_detectors_HEADERS
        PRIVATE detray::core_array detray::io_array
    )
endif()

# Shipped names use the manifest; external files reuse the same generator.
foreach(GENERATOR IN LISTS DETRAY_GENERATE_METADATA)
    detray_add_metadata_header(detray_detectors "${GENERATOR}")
endforeach()
