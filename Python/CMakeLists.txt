acts_apply_warning_flags(CXX ${ACTS_FAIL_ON_WARNINGS})

# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# List to collect built plugins
set(_plugins_built)

# Shared extension, installation, environment, and stub-generation helpers.
include(ActsPython)
acts_python_environment(${_python_dir} python)

# Add the bindings
add_subdirectory(Core)
add_subdirectory_if(Fatras ACTS_BUILD_FATRAS)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory(Plugins)
add_subdirectory(Utilities)

# Install all generated stubs in one pass
acts_python_install_stubs(
    ${_python_dir}/acts
    ${_acts_python_install_dir}
    ${ACTS_GENERATE_PYTHON_STUBS}
)
