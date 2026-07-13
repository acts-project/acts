# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir}/detray)

if(DETRAY_BUILD_TEST_UTILS)
    pybind11_add_module(DetrayExamplesPythonBindings src/ExamplesModuleEntry.cpp)

    target_link_libraries(
        DetrayExamplesPythonBindings
        PRIVATE
            detray::core_array
            detray::io
            detray::test_common
            detray::test_utils
    )

    # Place the compiled module inside the detray package directory.
    set_target_properties(
        DetrayExamplesPythonBindings
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/detray
    )
    install(
        TARGETS DetrayExamplesPythonBindings
        DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}
    )

    set(py_files examples.py)

    foreach(f ${py_files})
        set(_target ${_python_dir}/detray/${f})
        get_filename_component(_dir ${_target} DIRECTORY)
        file(MAKE_DIRECTORY ${_dir})

        file(
            CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/python/${f} ${_target}
            SYMBOLIC
        )

        get_filename_component(_rel ${f} DIRECTORY)

        install(
            FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/${f}
            DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}/${_rel}
        )
    endforeach()

    detray_pybind_generate_stubs(detray.examples DetrayExamplesPythonBindings)
endif()
