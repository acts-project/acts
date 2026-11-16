# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

if(NOT DETRAY_BUILD_PYTHON_BINDINGS)
    return()
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/ActsPython.cmake)
set(_python_dir "${CMAKE_BINARY_DIR}/python")
acts_python_environment(${_python_dir} python)
acts_python_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/package_init.py.in
    ${_python_dir}/detray/__init__.py
    ${DETRAY_PYTHON_INSTALL_DIR}
)

add_library(detray_python_common INTERFACE)
add_library(detray::python_common ALIAS detray_python_common)
target_include_directories(
    detray_python_common
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(
    detray_python_common
    INTERFACE detray::core_array detray::detectors pybind11::headers
)

add_subdirectory(Core)
add_subdirectory(Io)
add_subdirectory(Plugins)
add_subdirectory(Tests)

acts_python_install_stubs(
    ${_python_dir}/detray
    ${DETRAY_PYTHON_INSTALL_DIR}
    ${DETRAY_GENERATE_PYTHON_STUBS}
)
