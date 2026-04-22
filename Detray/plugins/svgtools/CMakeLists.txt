# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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

# Clean up.
unset(_detray_svgtools_public_headers)
