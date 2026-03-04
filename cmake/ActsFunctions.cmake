# This function creates a library with Acts namespace support
# It takes the library name without the "Acts" prefix and creates:
# 1. A direct CMake target named "Acts${library_name}"
# 2. An ALIAS target named "Acts::${library_name}"
function(acts_add_library library_name)
    set(options "")
    set(oneValueArgs ACTS_INCLUDE_FOLDER)
    set(multiValueArgs "")
    cmake_parse_arguments(
        PARSE_ARGV
        1
        lib_args
        "${options}"
        "${oneValueArgs}"
        "${multiValueArgs}"
    )

    set(full_target_name "Acts${library_name}")

    # Forward all arguments except the first (library_name) to add_library
    add_library(${full_target_name} ${lib_args_UNPARSED_ARGUMENTS})

    # Create alias target with Acts:: namespace
    add_library(Acts::${library_name} ALIAS ${full_target_name})

    set_target_properties(
        ${full_target_name}
        PROPERTIES EXPORT_NAME Acts::${library_name}
    )

    install(
        TARGETS ${full_target_name}
        EXPORT ${full_target_name}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    if(lib_args_ACTS_INCLUDE_FOLDER)
        install(
            DIRECTORY ${lib_args_ACTS_INCLUDE_FOLDER}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    endif()
endfunction()

# This function adds a helper target called ${target}_HEADERS which has
# generated source files that include each hedaer that was given one-by-one.
# The generated target links against the one given in `target` and therefore
# has access to all includes.
# The generated target is not included in the default build, and is only meant
# to help tools like clangd discover compiler flags.
function(acts_compile_headers target)
    set(options "")
    set(oneValueArgs "")
    set(multiValueArgs GLOB)
    cmake_parse_arguments(
        PARSE_ARGV
        0
        ARGS
        "${options}"
        "${oneValueArgs}"
        "${multiValueArgs}"
    )

    if(NOT ACTS_COMPILE_HEADERS)
        return()
    endif()

    if(NOT "${ARGS_GLOB}" STREQUAL "")
        if(NOT "${_headers}" STREQUAL "")
            message(SEND_ERROR "Cannot use HEADERS and GLOB at the same time")
            return()
        endif()

        file(
            GLOB_RECURSE _headers
            RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
            ${ARGS_GLOB}
        )
    endif()

    if("${_headers}" STREQUAL "")
        message(SEND_ERROR "No headers specified")
        return()
    endif()

    set(_sources "")

    foreach(_file ${_headers})
        if(IS_ABSOLUTE "${_file}")
            message(SEND_ERROR "Absolute paths are not allowed: ${_file}")
            continue()
        endif()

        set(_header_file "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")

        if(NOT EXISTS "${_header_file}")
            message(SEND_ERROR "File not found: ${_header_file}")
        endif()
        if(IS_DIRECTORY "${_header_file}")
            message(SEND_ERROR "Path is a directory: ${_header_file}")
        endif()

        get_filename_component(_header_file_name "${_file}" NAME)
        get_filename_component(_header_directory "${_file}" DIRECTORY)

        set(_temporary_dir "${CMAKE_CURRENT_BINARY_DIR}/${_header_directory}")

        file(MAKE_DIRECTORY "${_temporary_dir}")

        set(_temporary_path "${_temporary_dir}/${_header_file_name}.cpp")

        file(WRITE "${_temporary_path}" "#include \"${_header_file}\"")

        list(APPEND _sources "${_temporary_path}")
    endforeach()

    if(NOT TARGET Acts${target}_HEADERS)
        add_library(Acts${target}_HEADERS SHARED EXCLUDE_FROM_ALL ${_sources})
        target_link_libraries(Acts${target}_HEADERS PRIVATE Acts::${target})
    else()
        target_sources(Acts${target}_HEADERS PRIVATE ${_sources})
    endif()
endfunction()
