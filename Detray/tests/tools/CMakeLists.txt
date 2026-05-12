# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building 'detray::tools' component")

include(CMakeFindDependencyMacro)

find_dependency(Boost COMPONENTS program_options)

add_library(detray_tools INTERFACE)
add_library(detray::tools ALIAS detray_tools)

target_include_directories(
    detray_tools
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    detray_tools
    INTERFACE Boost::program_options vecmem::core detray::io detray::test_common
)

add_subdirectory(src)
