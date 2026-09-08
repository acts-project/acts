# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Let the user know what's happening.
message(STATUS "Building 'detray::svgtools' component")

# Set up the library.
file(
    GLOB _detray_svgtools_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/plugins/svgtools/*.hpp"
)

detray_add_library( detray_svgtools svgtools
   ${_detray_svg_public_headers}
)
target_link_libraries(
    detray_svgtools
    INTERFACE detray::core detray::io actsvg::core actsvg::meta
)

# Compile all headers (including the ones under detail/) one-by-one.
detray_compile_headers(
    svgtools
    detray::svgtools
    "include/detray/plugins/svgtools/*.hpp"
    "include/detray/plugins/svgtools/**/*.hpp"
)
if(TARGET detray_svgtools_HEADERS)
    # The svgtools headers need an algebra plugin to compile.
    target_link_libraries(detray_svgtools_HEADERS PRIVATE detray::core_array)
endif()

# Clean up.
unset(_detray_svgtools_public_headers)
