# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# List to collect built plugins
set(_plugins_built)

# Generate stubs for a user-facing Python module using pybind11-stubgen.
# _module:       dotted Python module name (e.g. acts, acts.examples, acts.json)
# _cmake_target: CMake target to attach POST_BUILD to (must be in the calling directory)
function(acts_pybind_generate_stubs _module _cmake_target)
    if(NOT ACTS_GENERATE_PYTHON_STUBS)
        return()
    endif()
    add_custom_command(
        TARGET ${_cmake_target}
        POST_BUILD
        COMMAND
            ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
            ${uv_exe} tool run --python ${Python_EXECUTABLE} --with numpy
            pybind11-stubgen ${_module} --exit-code --output-dir
            ${CMAKE_BINARY_DIR}/python --ignore-invalid-expressions ".*"
            --ignore-invalid-identifiers "GeoModelTree::FpvConstLink"
            --ignore-unresolved-names ".*"
        COMMENT "Generating stubs for ${_module}"
        VERBATIM
    )
endfunction()

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
