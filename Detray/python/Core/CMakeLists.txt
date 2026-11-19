# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

set(_python_dir "${CMAKE_BINARY_DIR}/python")
set(_pkg_dir ${_python_dir}/detray/core)
file(MAKE_DIRECTORY ${_pkg_dir})

acts_add_python_module(
    DetrayPythonBindings
    SOURCES src/CoreModuleEntry.cpp
    OUTPUT_DIRECTORY ${_pkg_dir}
    INSTALL_DIRECTORY ${DETRAY_PYTHON_INSTALL_DIR}/core
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)

target_link_libraries(
    DetrayPythonBindings
    PRIVATE detray::core_array detray::io detray::detectors vecmem::core
)

# Generate the package __init__ that re-exports the compiled module.
set(BINDINGS_IMPORTS "")
set(BINDINGS_MODULE DetrayPythonBindings)
acts_python_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/../module_init.py.in
    ${_pkg_dir}/__init__.py
    ${DETRAY_PYTHON_INSTALL_DIR}/core
)

acts_pybind_generate_stubs(
    detray.core
    DetrayPythonBindings
    PACKAGE
    ENABLED ${DETRAY_GENERATE_PYTHON_STUBS}
    PYTHON_ROOT ${_python_dir}
    PYTHON_EXECUTABLE ${Python_EXECUTABLE}
)
