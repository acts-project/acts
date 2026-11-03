# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

set(_python_dir "${CMAKE_BINARY_DIR}/python")
set(_pkg_dir ${_python_dir}/detray/tests)

if(NOT TARGET detray::svgtools)
    message(
        STATUS
        "Skipping 'detray.tests' python module: 'detray::svgtools' not available"
    )
    return()
endif()

file(MAKE_DIRECTORY ${_pkg_dir})

pybind11_add_module(DetrayTestsPythonBindings src/TestsModuleEntry.cpp)

target_link_libraries(
    DetrayTestsPythonBindings
    PRIVATE
        detray::core_array
        detray::io
        detray::test_common
        detray::svgtools
        GTest::gtest
        vecmem::core
        covfie::core
)

# Place the compiled module inside the `detray.tests` package.
set_target_properties(
    DetrayTestsPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_pkg_dir}
)
install(
    TARGETS DetrayTestsPythonBindings
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/tests
)

# Generate the package __init__ that re-exports the compiled module.
set(BINDINGS_MODULE DetrayTestsPythonBindings)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/module_init.py.in
    ${_pkg_dir}/__init__.py
    @ONLY
)
install(
    FILES ${_pkg_dir}/__init__.py
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/tests
)

acts_pybind_generate_stubs(detray.tests DetrayTestsPythonBindings)
