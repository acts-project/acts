# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)
set(_python_install_dir "python/acts")
set(_plugins_built)

# Add the bindings
add_subdirectory(Utilities)
add_subdirectory(Core)
add_subdirectory(Plugins)

# Finish up the __init__ file with the correct modules built
# This will make them available through import * from acts
string(REPLACE ";" "," _plugins_built "${_plugins_built}")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/python/__init__.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/acts/__init__.py
    @ONLY
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/acts/__init__.py
    DESTINATION ${_python_install_dir}
)
