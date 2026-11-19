include_guard(GLOBAL)
include(${CMAKE_CURRENT_LIST_DIR}/ActsPybindStubs.cmake)

# Shared by the ACTS and standalone detray entry points. Callers own package
# layout and feature policy; this module owns extension/install mechanics.
function(acts_add_python_module target)
    cmake_parse_arguments(
        PARSE_ARGV
        1
        ARG
        ""
        "OUTPUT_DIRECTORY;INSTALL_DIRECTORY;LIBRARY_DIRECTORY"
        "SOURCES"
    )
    if(ARG_UNPARSED_ARGUMENTS OR ARG_KEYWORDS_MISSING_VALUES)
        message(
            FATAL_ERROR
            "Invalid arguments to acts_add_python_module(${target})"
        )
    endif()
    foreach(
        _required
        OUTPUT_DIRECTORY
        INSTALL_DIRECTORY
        LIBRARY_DIRECTORY
        SOURCES
    )
        if(NOT DEFINED ARG_${_required})
            message(FATAL_ERROR "acts_add_python_module requires ${_required}")
        endif()
    endforeach()
    pybind11_add_module(${target} ${ARG_SOURCES})
    if(APPLE)
        set(_origin "@loader_path")
    else()
        set(_origin "$ORIGIN")
    endif()
    get_filename_component(
        _package
        "${ARG_INSTALL_DIRECTORY}"
        ABSOLUTE
        BASE_DIR "${CMAKE_INSTALL_PREFIX}"
    )
    get_filename_component(
        _library
        "${ARG_LIBRARY_DIRECTORY}"
        ABSOLUTE
        BASE_DIR "${CMAKE_INSTALL_PREFIX}"
    )
    file(RELATIVE_PATH _relative "${_package}" "${_library}")
    set_target_properties(
        ${target}
        PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${ARG_OUTPUT_DIRECTORY}"
            INSTALL_RPATH "${_origin}/${_relative}"
    )
    install(TARGETS ${target} DESTINATION "${ARG_INSTALL_DIRECTORY}")
endfunction()

# Copy/configure a Python source once into the build tree and install that same
# file. Variables used by templates are inherited from the caller.
function(acts_python_file source output install_directory)
    if(source MATCHES "\\.in$")
        configure_file("${source}" "${output}" @ONLY)
    else()
        configure_file("${source}" "${output}" COPYONLY)
    endif()
    install(FILES "${output}" DESTINATION "${install_directory}")
endfunction()

# The two entry points may share an output tree. Only one setup script is needed.
function(acts_python_environment output_directory install_directory)
    get_property(_outputs GLOBAL PROPERTY ACTS_PYTHON_ENVIRONMENTS)
    if(output_directory IN_LIST _outputs)
        return()
    endif()
    set_property(
        GLOBAL
        APPEND
        PROPERTY ACTS_PYTHON_ENVIRONMENTS "${output_directory}"
    )
    acts_python_file(
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/python_setup.sh.in"
        "${output_directory}/setup.sh"
        "${install_directory}"
    )
endfunction()
