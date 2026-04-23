# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# List to collect built plugins
set(_plugins_built)

function(acts_pybind_generate_stubs _module _install_dest)
    if(NOT ACTS_GENERATE_PYTHON_STUBS)
        return()
    endif()
    string(REPLACE "." "/" _module_path "${_module}")
    add_custom_command(
        TARGET ActsPythonBindings POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
            ${uv_exe} tool run --python ${Python_EXECUTABLE} pybind11-stubgen ${_module}
            --output-dir ${CMAKE_BINARY_DIR}/python
            --ignore-invalid-expressions ".*"
            --ignore-invalid-identifiers "GeoModelTree::FpvConstLink"
            --ignore-unresolved-names "_abc\\._abc_data"
        COMMENT "Generating stubs for ${_module}"
        VERBATIM
    )
    install(
        DIRECTORY ${CMAKE_BINARY_DIR}/python/${_module_path}/
        DESTINATION ${_install_dest}
        FILES_MATCHING PATTERN "*.pyi"
    )
endfunction()

# Add the bindings
add_subdirectory(Core)
add_subdirectory(Fatras)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory(Plugins)
add_subdirectory(Utilities)

