include_guard(GLOBAL)

if(NOT "${CMAKE_CURRENT_LIST_DIR}" IN_LIST CMAKE_MODULE_PATH)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
endif()

# Generate stubs for a user-facing Python module using pybind11-stubgen.
# _module:       dotted Python module name (e.g. acts, acts.examples, detray.core)
# _cmake_target: primary CMake target whose .so must exist before stub generation
# [extra_deps]:  additional CMake targets that must be built first (e.g. ActsPythonBindings)
function(acts_pybind_generate_stubs _module _cmake_target)
    if(NOT ACTS_GENERATE_PYTHON_STUBS)
        return()
    endif()
    include(ActsEnsureUv)
    set(_extra_deps ${ARGN})
    string(REPLACE "." "/" _module_path ${_module})
    set(_stub_file "${CMAKE_BINARY_DIR}/python/${_module_path}/__init__.pyi")
    add_custom_command(
        OUTPUT ${_stub_file}
        COMMAND
            ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_BINARY_DIR}/python"
            ${uv_exe} tool run --python ${Python_EXECUTABLE} --with numpy
            pybind11-stubgen ${_module} --exit-code --output-dir
            ${CMAKE_BINARY_DIR}/python --ignore-invalid-expressions ".*"
            --ignore-invalid-identifiers "GeoModelTree::FpvConstLink"
            --ignore-unresolved-names ".*"
        DEPENDS ${_cmake_target} ${_extra_deps}
        COMMENT "Generating stubs for ${_module}"
        VERBATIM
    )
    add_custom_target(${_cmake_target}Stubs ALL DEPENDS ${_stub_file})
endfunction()
