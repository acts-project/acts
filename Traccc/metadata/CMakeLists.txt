# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Detector metadata headers that traccc needs but detray does not ship. They are
# generated from external detray metadata configurations and installed with
# traccc.
# traccc::itk_detector includes a fixed header and uses detray::itk_metadata.
# Reject selectors that would generate a different name at configure time.
if(
    NOT TRACCC_ITK_METADATA_CONFIGURATION
        MATCHES
        "(^|[/\\\\])itk(_detector)?\\.py:[A-Za-z_][A-Za-z0-9_]*$"
)
    message(
        FATAL_ERROR
        "TRACCC_ITK_METADATA_CONFIGURATION must select an external ITk configuration "
        "(/path/to/itk_detector.py:configure or /path/to/itk.py:configure)"
    )
endif()

traccc_add_library(traccc_detray_metadata detray_metadata TYPE INTERFACE)
target_link_libraries(traccc_detray_metadata INTERFACE detray::detectors)
target_include_directories(
    traccc_detray_metadata
    INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# detray/detectors/itk_metadata.hpp, used by traccc::itk_detector.
detray_add_metadata_header(
    traccc_detray_metadata
    "${TRACCC_ITK_METADATA_CONFIGURATION}"
)
