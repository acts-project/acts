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

# Generate the metadata headers. Only hand-written headers (telescope) live in
# the source tree; the shipped detector configurations in
# codegen/detray_metadata/configurations are always generated.
set(_detray_shipped_metadata default odd toy wire_chamber)
foreach(_name IN LISTS _detray_shipped_metadata)
    detray_add_metadata_header(detray_detectors "${_name}")
endforeach()

# DETRAY_GENERATE_METADATA only adds external configurations.
foreach(_configuration IN LISTS DETRAY_GENERATE_METADATA)
    string(REGEX REPLACE "_metadata$" "" _name "${_configuration}")
    if(_name IN_LIST _detray_shipped_metadata)
        message(
            DEPRECATION
            "DETRAY_GENERATE_METADATA: '${_configuration}' is always generated, "
            "remove it from the list"
        )
        continue()
    endif()
    if(_name STREQUAL "itk")
        message(
            DEPRECATION
            "DETRAY_GENERATE_METADATA: the ITk metadata is no longer shipped by "
            "detray; traccc generates it (see TRACCC_ITK_METADATA_CONFIGURATION)"
        )
        continue()
    endif()
    detray_add_metadata_header(detray_detectors "${_configuration}")
endforeach()

unset(_detray_shipped_metadata)
