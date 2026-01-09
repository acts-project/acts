# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# If we are in a scikit-build-core build, install to SKBUILD_PLATLIB_DIR
# Otherwise, we install to CMAKE_INSTALL_PREFIX/python
if(ACTS_BUILD_WITH_SKBUILD)
    set(ACTS_PYTHON_INSTALL_PREFIX ${SKBUILD_PLATLIB_DIR})
else()
    set(ACTS_PYTHON_INSTALL_PREFIX "python")
endif()

set(_python_install_dir "${ACTS_PYTHON_INSTALL_PREFIX}/acts")

set(_plugins_built)

# Add the bindings
add_subdirectory(Core)
add_subdirectory(Examples)
add_subdirectory(Plugins)
add_subdirectory(Utilities)
