# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# CMake include(s).
include(CMakeParseArguments)

# Helper function to create and set the `--keep` flag on CUDA targets.
function(detray_add_cuda_artifact_dir_to_target target)
    set(cuda_keep_dir
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/cuda_artifacts/${target}/"
    )
    add_custom_target(
        ${target}_cuda_artifact_mkdir
        COMMAND ${CMAKE_COMMAND} -E make_directory ${cuda_keep_dir}
    )
    add_dependencies(${target} ${target}_cuda_artifact_mkdir)
    target_compile_options(
        ${target}
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--keep --keep-dir ${cuda_keep_dir}>
    )
endfunction(detray_add_cuda_artifact_dir_to_target)

# Helper function for setting up the detray libraries.
#
# Usage: detray_add_library( detray_core core "header1.hpp"... )
#
function(detray_add_library fullname basename)
    # Create the library.
    add_library(${fullname} INTERFACE ${ARG_UNPARSED_ARGUMENTS})

    # Set up how clients should find its headers.
    target_include_directories(
        ${fullname}
        INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    # Make sure that the library is available as "detray::${basename}" in every
    # situation.
    set_target_properties(${fullname} PROPERTIES EXPORT_NAME ${basename})
    add_library(detray::${basename} ALIAS ${fullname})

    # Set up the installation of the library and its headers.
    install(
        TARGETS ${fullname}
        EXPORT detray-exports
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        OPTIONAL
    )
endfunction(detray_add_library)

# Helper function testing the detray public headers.
#
# It can be used to test that public headers would include everything
# that they need to work, and that the CMake library targets would take
# care of declaring all of their dependencies correctly for the public
# headers to work.
#
# Usage: detray_test_public_headers( detray_core
#                                    include/header1.hpp ... )
#
function(detray_test_public_headers library)
    # If testing is not turned on, don't do anything.
    if((NOT BUILD_TESTING) OR (NOT DETRAY_BUILD_TESTING))
        return()
    endif()

    # All arguments are treated as header file names.
    foreach(_headerName ${ARGN})
        # Make the header filename into a "string".
        string(REPLACE "/" "_" _headerNormName "${_headerName}")
        string(REPLACE "." "_" _headerNormName "${_headerNormName}")

        # Write a small source file that would test that the public
        # header can be used as-is.
        set(_testFileName
            "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/test_${library}_${_headerNormName}.cpp"
        )
        if(NOT EXISTS "${_testFileName}")
            file(
                WRITE "${_testFileName}"
                "#include \"${_headerName}\"\n"
                "int main() { return 0; }"
            )
        endif()

        # Set up an executable that would build it. But hide it, don't put it
        # into ${CMAKE_BINARY_DIR}/bin.
        add_executable("test_${library}_${_headerNormName}" "${_testFileName}")
        target_link_libraries(
            "test_${library}_${_headerNormName}"
            PRIVATE ${library}
        )
        set_target_properties(
            "test_${library}_${_headerNormName}"
            PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY
                    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}"
        )
    endforeach()
endfunction(detray_test_public_headers)

# Helper function for setting up the detray executables.
#
# The detray executables are *not* installed with the project, as they are only
# used for testing / benchmarking the code. Clients of detray do not need them.
#
# Usage: detray_add_executable( foo bar.cpp
#                               LINK_LIBRARIES detray::core )
#
function(detray_add_executable name)
    # Parse the function's options.
    cmake_parse_arguments(ARG "" "" "LINK_LIBRARIES" ${ARGN})

    # Create the executable.
    set(exe_name "detray_${name}")
    add_executable(${exe_name} ${ARG_UNPARSED_ARGUMENTS})
    if(ARG_LINK_LIBRARIES)
        target_link_libraries(${exe_name} PRIVATE ${ARG_LINK_LIBRARIES})
    endif()

    detray_add_cuda_artifact_dir_to_target(${exe_name})
endfunction(detray_add_executable)

# Helper function for setting up the detray tests.
#
# Usage: detray_add_test( core source1.cpp source2.cpp
#                         LINK_LIBRARIES detray::core )
#
function(detray_add_test name)
    # Parse the function's options.
    cmake_parse_arguments(ARG "" "" "LINK_LIBRARIES" ${ARGN})

    # Create the test executable.
    set(test_exe_name "detray_${name}")
    add_executable(${test_exe_name} ${ARG_UNPARSED_ARGUMENTS})
    if(ARG_LINK_LIBRARIES)
        target_link_libraries(${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES})
    endif()

    # Run tests with sanitizers
    if(DETRAY_ENABLE_SANITIZER)
        # Common flags
        set(SANITIZER_FLAGS
            "-fsanitize=address,undefined,pointer-compare,pointer-subtract,float-divide-by-zero"
        )

        # Extra flags for clang
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(SANITIZER_FLAGS
                "${SANITIZER_FLAGS},leak,integer,nullability,implicit-conversion,local-bounds"
            )
        endif()

        target_compile_options(${test_exe_name} PUBLIC ${SANITIZER_FLAGS})
        target_link_options(${test_exe_name} PUBLIC ${SANITIZER_FLAGS})

        # Clean up
        unset(SANITIZER_FLAGS)
    endif()

    # Run the executable as the test.
    add_test(NAME ${test_exe_name} COMMAND ${test_exe_name})

    # Set all properties for the test.
    set_tests_properties(
        ${test_exe_name}
        PROPERTIES ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data/
    )

    detray_add_cuda_artifact_dir_to_target(${test_exe_name})
endfunction(detray_add_test)

# Helper function to set up a unit test
#
# Usage: See detray_add_test
#
function(detray_add_unit_test name)
    detray_add_test(unit_test_${name} ${ARGN})
endfunction(detray_add_unit_test)

# Helper function to set up an integration test
#
# Usage: See detray_add_test
#
function(detray_add_integration_test name)
    detray_add_test(integration_test_${name} ${ARGN})
endfunction(detray_add_integration_test)

# Helper function for setting up the detray tutorials.
#
# Usage: detray_add_tutorial( core source1.cpp source2.cpp
#                             LINK_LIBRARIES detray::core )
#
function(detray_add_tutorial name)
    # Parse the function's options.
    cmake_parse_arguments(ARG "" "" "LINK_LIBRARIES" ${ARGN})

    # Create the tutorial executable.
    set(tutorial_exe_name "detray_tutorial_${name}")
    add_executable(${tutorial_exe_name} ${ARG_UNPARSED_ARGUMENTS})
    if(ARG_LINK_LIBRARIES)
        target_link_libraries(
            ${tutorial_exe_name}
            PRIVATE ${ARG_LINK_LIBRARIES}
        )
    endif()

    detray_add_cuda_artifact_dir_to_target(${tutorial_exe_name})
endfunction(detray_add_tutorial)

# Helper function for adding individual flags to "flag variables".
#
# Usage: detray_add_flag( CMAKE_CXX_FLAGS "-Wall" )
#
function(detray_add_flag name value)
    # Escape special characters in the value:
    set(matchedValue "${value}")
    foreach(
        c
        "*"
        "."
        "^"
        "$"
        "+"
        "?"
    )
        string(REPLACE "${c}" "\\${c}" matchedValue "${matchedValue}")
    endforeach()

    # Check if the variable already has this value in it:
    if("${${name}}" MATCHES "${matchedValue}")
        return()
    endif()

    # If not, then let's add it now:
    set(${name} "${${name}} ${value}" PARENT_SCOPE)
endfunction(detray_add_flag)

# Generate a single detray sympy codegen header via acts_code_generation and
# install it alongside the other detray headers.
#
# Usage: detray_add_codegen_header(
#            TARGET detray_core
#            CODEGEN_DIR ${detray_codegen_dir}
#            ACTS_CODEGEN_PKG ${acts_codegen_pkg}
#            SCRIPT gen_full_jacobian.py
#            OUTPUT detray/propagator/actors/codegen/full_jacobian.hpp
#        )
#
function(detray_add_codegen_header)
    set(oneValueArgs TARGET CODEGEN_DIR ACTS_CODEGEN_PKG SCRIPT OUTPUT)
    cmake_parse_arguments(ARG "" "${oneValueArgs}" "" ${ARGN})

    if(
        NOT ARG_TARGET
        OR NOT ARG_SCRIPT
        OR NOT ARG_OUTPUT
        OR NOT ARG_CODEGEN_DIR
        OR NOT ARG_ACTS_CODEGEN_PKG
    )
        message(
            FATAL_ERROR
            "detray_add_codegen_header: TARGET, CODEGEN_DIR, ACTS_CODEGEN_PKG, SCRIPT and OUTPUT are required"
        )
    endif()

    acts_code_generation(
        ADD_TO_TARGET ${ARG_TARGET}
        PYTHON ${ARG_CODEGEN_DIR}/${ARG_SCRIPT}
        WITH_REQUIREMENTS ${ARG_ACTS_CODEGEN_PKG}/requirements.txt
        WITH ${ARG_ACTS_CODEGEN_PKG}
        ISOLATED
        OUTPUT ${ARG_OUTPUT}
        RESULT_INCLUDE_DIR _gen_root
    )

    get_filename_component(_output_subdir ${ARG_OUTPUT} DIRECTORY)
    install(
        FILES ${_gen_root}/${ARG_OUTPUT}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${_output_subdir}
    )
endfunction(detray_add_codegen_header)
