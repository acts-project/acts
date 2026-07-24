# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# List to collect built plugins
set(_plugins_built)

# Provides acts_pybind_generate_stubs() (shared with detray).
include(ActsPybindStubs)

# Add the bindings
add_subdirectory(Core)
add_subdirectory(Fatras)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory(Plugins)
add_subdirectory(Utilities)

# Install all generated stubs in one pass
if(ACTS_GENERATE_PYTHON_STUBS)
    install(
        DIRECTORY ${CMAKE_BINARY_DIR}/python/acts/
        DESTINATION ${_acts_python_install_dir}
        FILES_MATCHING
        PATTERN "*.pyi"
    )
endif()
