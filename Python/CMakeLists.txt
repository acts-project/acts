set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)
set(_python_install_dir "python/acts")

# Add the utilities as ActsPythonUtilties
add_subdirectory(Utilities)

# This creates the ActsPythonBindings library
add_subdirectory(Module)

# Add the core module
add_subdirectory(Core)