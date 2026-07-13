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

pybind11_add_module(DetrayPythonBindings src/CoreModuleEntry.cpp)

target_link_libraries(
    DetrayPythonBindings
    PRIVATE
        detray::core_array
        detray::io
        detray::detectors
        vecmem::core
        covfie::core
)

# Place the compiled module inside the `detray.core` package.
set_target_properties(
    DetrayPythonBindings
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_pkg_dir}
)
install(
    TARGETS DetrayPythonBindings
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/core
)

set(py_files __init__.py)

foreach(f ${py_files})
    set(_target ${_pkg_dir}/${f})
    get_filename_component(_dir ${_target} DIRECTORY)
    file(MAKE_DIRECTORY ${_dir})

    file(
        CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
        SYMBOLIC
    )

    get_filename_component(_rel ${f} DIRECTORY)

    install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
        DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/core/${_rel}
    )
endforeach()

detray_pybind_generate_stubs(detray.core DetrayPythonBindings)
