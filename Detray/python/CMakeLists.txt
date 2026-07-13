# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the python stall directories
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/detray)

# Set up the source script in the python directory
configure_file(setup.sh.in ${_python_dir}/setup.sh @ONLY)
install(FILES ${_python_dir}/setup.sh DESTINATION python)

# Set up the init file for the detray python module
set(_target ${_python_dir}/detray/__init__.py)
file(
    CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/detray/__init__.py ${_target}
    SYMBOLIC
)

# Allow access to the python detectors package from the main __init__.py
file(
    CREATE_LINK
        ${CMAKE_CURRENT_SOURCE_DIR}/detray/detectors
        ${_python_dir}/detray/detectors
    SYMBOLIC
)

install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/detray/__init__.py
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}
)

if(NOT COMMAND acts_pybind_generate_stubs)
    include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/ActsPybindStubs.cmake)
endif()
set(ACTS_GENERATE_PYTHON_STUBS ${DETRAY_GENERATE_PYTHON_STUBS})

if(DETRAY_BUILD_PYTHON_BINDINGS)
    add_subdirectory(Core)
    add_subdirectory(Examples)
endif()

# Install the generated stub files alongside the package.
if(DETRAY_GENERATE_PYTHON_STUBS)
    install(
        DIRECTORY ${_python_dir}/detray/
        DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}
        FILES_MATCHING
        PATTERN "*.pyi"
    )
endif()
