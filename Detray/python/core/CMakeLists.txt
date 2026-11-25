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

pybind11_add_module(detray_core_python_bindings src/core_module_entry.cpp)

target_link_libraries(
    detray_core_python_bindings
    PRIVATE detray::core_array detray::detectors vecmem::core
)

# Place the compiled module inside the `detray.core` package.
set_target_properties(
    detray_core_python_bindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_pkg_dir}
)
install(
    TARGETS detray_core_python_bindings
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/core
)

# Generate the package __init__ that re-exports the compiled module.
set(BINDINGS_MODULE detray_core_python_bindings)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/../module_init.py.in
    ${_pkg_dir}/__init__.py
    @ONLY
)
install(
    FILES ${_pkg_dir}/__init__.py
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/core
)

acts_pybind_generate_stubs(detray.core detray_core_python_bindings)
