# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

set(_python_dir "${CMAKE_BINARY_DIR}/python")
set(_pkg_dir ${_python_dir}/detray/tests)

if(
    NOT TARGET detray::svgtools
    OR NOT TARGET detray::test_common
    OR NOT TARGET detray::test_framework
    OR NOT TARGET GTest::gtest
)
    message(
        STATUS
        "Skipping 'detray.tests' python module: validation dependencies not available"
    )
    return()
endif()

file(MAKE_DIRECTORY ${_pkg_dir})

acts_add_python_module(
    DetrayTestsPythonBindings
    SOURCES src/TestsModuleEntry.cpp
    OUTPUT_DIRECTORY ${_pkg_dir}
    INSTALL_DIRECTORY ${DETRAY_PYTHON_INSTALL_DIR}/tests
    LIBRARY_DIRECTORY ${CMAKE_INSTALL_LIBDIR}
)

target_link_libraries(
    DetrayTestsPythonBindings
    PRIVATE
        detray::python_common
        detray::core_array
        detray::io
        detray::test_common
        detray::test_framework
        detray::svgtools
        GTest::gtest
        vecmem::core
        covfie::core
)

# Validation is a public Python utility, so a locally built shared GoogleTest
# is a runtime dependency of the installed module, not only a build dependency.
get_target_property(_gtest_target GTest::gtest ALIASED_TARGET)
if(NOT _gtest_target)
    set(_gtest_target GTest::gtest)
endif()
get_target_property(_gtest_imported ${_gtest_target} IMPORTED)
get_target_property(_gtest_type ${_gtest_target} TYPE)
if(NOT _gtest_imported AND _gtest_type STREQUAL "SHARED_LIBRARY")
    install(
        TARGETS ${_gtest_target}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

# Generate the package __init__ that re-exports the compiled module.
set(BINDINGS_IMPORTS "from .. import core as _core\n\n")
set(BINDINGS_MODULE DetrayTestsPythonBindings)
acts_python_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/../module_init.py.in
    ${_pkg_dir}/__init__.py
    ${DETRAY_PYTHON_INSTALL_DIR}/tests
)

acts_pybind_generate_stubs(
    detray.tests
    DetrayTestsPythonBindings
    PACKAGE
    ENABLED ${DETRAY_GENERATE_PYTHON_STUBS}
    PYTHON_ROOT ${_python_dir}
    PYTHON_EXECUTABLE ${Python_EXECUTABLE}
    DEPENDS DetrayPythonBindings
)
