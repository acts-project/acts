include_guard(GLOBAL)

# Generate stubs with explicit package policy, independent of the entry point.
function(acts_pybind_generate_stubs _module _cmake_target)
    cmake_parse_arguments(
        PARSE_ARGV
        2
        ARG
        "PACKAGE"
        "ENABLED;PYTHON_ROOT;PYTHON_EXECUTABLE"
        "DEPENDS;IGNORE_INVALID_IDENTIFIERS"
    )
    if(ARG_UNPARSED_ARGUMENTS OR ARG_KEYWORDS_MISSING_VALUES)
        message(FATAL_ERROR "Invalid arguments to acts_pybind_generate_stubs")
    endif()
    if(NOT DEFINED ARG_ENABLED)
        message(
            FATAL_ERROR
            "Stub generation requires an explicit ENABLED setting"
        )
    endif()
    if(NOT ARG_ENABLED)
        return()
    endif()
    if(NOT ARG_PYTHON_ROOT OR NOT ARG_PYTHON_EXECUTABLE)
        message(
            FATAL_ERROR
            "Stub generation requires PYTHON_ROOT and PYTHON_EXECUTABLE"
        )
    endif()
    include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/ActsEnsureUv.cmake)
    set(_ignore_args)
    foreach(_identifier IN LISTS ARG_IGNORE_INVALID_IDENTIFIERS)
        list(APPEND _ignore_args --ignore-invalid-identifiers "${_identifier}")
    endforeach()
    string(REPLACE "." "/" _module_path ${_module})
    if(ARG_PACKAGE)
        set(_module_source "${ARG_PYTHON_ROOT}/${_module_path}/__init__.py")
        set(_stub_file "${_module_source}i")
    else()
        set(_module_source "${ARG_PYTHON_ROOT}/${_module_path}.py")
        set(_stub_file "${_module_source}i")
    endif()
    add_custom_command(
        OUTPUT ${_stub_file}
        COMMAND
            ${CMAKE_COMMAND} -E env "PYTHONPATH=${ARG_PYTHON_ROOT}" ${uv_exe}
            tool run --python ${ARG_PYTHON_EXECUTABLE} --with numpy
            pybind11-stubgen ${_module} --exit-code --output-dir
            ${ARG_PYTHON_ROOT} --ignore-invalid-expressions ".*" ${_ignore_args}
            --ignore-unresolved-names ".*"
        DEPENDS ${_cmake_target} ${ARG_DEPENDS} "${_module_source}"
        COMMENT "Generating stubs for ${_module}"
        VERBATIM
    )
    add_custom_target(${_cmake_target}Stubs ALL DEPENDS ${_stub_file})
endfunction()

# Stubgen can emit additional stubs for compiled submodules. Install the complete
# package's stub tree, using the same rules in both entry points.
function(acts_python_install_stubs source_directory install_directory enabled)
    if(enabled)
        install(
            DIRECTORY "${source_directory}/"
            DESTINATION "${install_directory}"
            FILES_MATCHING
            PATTERN "*.pyi"
        )
    endif()
endfunction()
