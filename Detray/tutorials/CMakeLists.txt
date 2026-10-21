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

# Set up a tutorial library, which showcases how to use detray components.
file(
    GLOB _detray_tutorial_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/tutorials/*.hpp"
)

detray_add_library( detray_tutorial tutorial ${_detray_tutorial_headers})
target_link_libraries(
    detray_tutorial
    INTERFACE vecmem::core detray::core detray::detectors
)
unset(_detray_tutorial_headers)

# Include all of the code-holding sub-directories.
add_subdirectory(src)
