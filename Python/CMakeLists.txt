# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)
set(_python_install_dir "python/acts")
set(_plugins_built)

# Add the bindings
add_subdirectory(Core)
add_subdirectory(Examples)
add_subdirectory(Plugins)
add_subdirectory(Utilities)
