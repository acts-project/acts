# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/detray/python/__init__.py
    DESTINATION ${DETRAY_PYTHON_INSTALL_DIR}
)
