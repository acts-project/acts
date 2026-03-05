include_guard(GLOBAL)

if(NOT ACTS_USE_SYSTEM_LIBS)
    message(STATUS "Configuring codegen: preparing uv")

    find_program(uv_exe uv)

    set(_uv_version "0.7.19")
    set(_base_url
        "https://github.com/astral-sh/uv/releases/download/${_uv_version}"
    )

    if(uv_exe STREQUAL "uv_exe-NOTFOUND")
        message(STATUS "uv not found, installing it")

        if(NOT APPLE AND NOT UNIX)
            message(FATAL_ERROR "Unsupported platform: ${CMAKE_SYSTEM_NAME}")
        endif()

        if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
            if(APPLE)
                set(UV_NAME "${_base_url}/uv-x86_64-apple-darwin.tar.gz")
                set(UV_HASH
                    "SHA256=698d24883fd441960fb4bc153b7030b89517a295502017ff3fdbba2fb0a0aa67"
                )
            elseif(UNIX)
                set(UV_URL "${_base_url}/uv-x86_64-unknown-linux-musl.tar.gz")
                set(UV_HASH
                    "SHA256=6236ed00a7442ab2c0f56f807d5a3331f3fb5c7640a357482fbc8492682641b2"
                )
            endif()
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "(arm)|(ARM)|(aarch64)")
            if(APPLE)
                set(UV_URL "${_base_url}/uv-aarch64-apple-darwin.tar.gz")
                set(UV_HASH
                    "SHA256=698d24883fd441960fb4bc153b7030b89517a295502017ff3fdbba2fb0a0aa67"
                )
            elseif(UNIX)
                set(UV_URL "${_base_url}/uv-aarch64-unknown-linux-musl.tar.gz")
                set(UV_HASH
                    "SHA256=e83c7c6d86c8e7456078c736a72550ce20222df8083f9317fc58cd49422ce5eb"
                )
            endif()
        else()
            message(
                FATAL_ERROR
                "Unsupported architecture: ${CMAKE_SYSTEM_PROCESSOR}"
            )
        endif()

        message(STATUS "Downloading uv from ${UV_URL}")
        set(UV_DIR "${CMAKE_BINARY_DIR}/uv")
        file(DOWNLOAD ${UV_URL} ${UV_DIR}/uv.tar.gz EXPECTED_HASH ${UV_HASH})

        file(ARCHIVE_EXTRACT INPUT ${UV_DIR}/uv.tar.gz DESTINATION ${UV_DIR})

        file(REMOVE ${UV_DIR}/uv.tar.gz)

        file(GLOB uv_extracted ${UV_DIR}/uv*)
        message(STATUS "Extracted uv: ${uv_extracted}")

        find_program(uv_exe uv PATHS ${uv_extracted} REQUIRED NO_DEFAULT_PATH)
    endif()

    message(STATUS "Found uv: ${uv_exe}")

    execute_process(
        COMMAND ${uv_exe} --version
        OUTPUT_VARIABLE uv_version
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "uv version: ${uv_version}")
else()
    message(
        STATUS
        "Configuring codegen in offline mode: preparing virtual environment"
    )

    find_package(Python REQUIRED COMPONENTS Interpreter)

    # The idea of the following code is to create a "nested" Python
    # environment; we grab the source of the packages in the current env
    # whether that is a virtual environment, a system environment, or a Spack
    # environment, and copy that into a newly created virtual environment.
    # This strategy comes from https://stackoverflow.com/a/75545634
    # First, we grab the Python package directory for the outside environment.
    execute_process(
        COMMAND
            ${Python_EXECUTABLE} -c
            "import sysconfig; print(sysconfig.get_paths()['purelib'])"
        OUTPUT_VARIABLE _python_package_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Then we create a new virtual env using the venv package which is built
    # into Python these days.
    execute_process(
        COMMAND ${Python_EXECUTABLE} -m venv ${CMAKE_BINARY_DIR}/codegen_venv
    )
    # Now, we get the package directory for the newly created virtual
    # environment.
    execute_process(
        COMMAND
            ${CMAKE_BINARY_DIR}/codegen_venv/bin/python -c
            "import sysconfig; print(sysconfig.get_paths()['purelib'])"
        OUTPUT_VARIABLE _python_nested_package_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Finally, we write the path found in the outside virtual env into the
    # new virtual env as described in the StackOverflow answer.
    file(
        WRITE
        "${_python_nested_package_dir}/_base_packages.pth"
        ${_python_package_dir}
    )

    message(
        STATUS
        "Virtual environment based on ${_python_package_dir} created in ${CMAKE_BINARY_DIR}/codegen_venv/"
    )
endif()

function(acts_code_generation)
    set(options ISOLATED)
    set(oneValueArgs ADD_TO_TARGET PYTHON PYTHON_VERSION OUTPUT)
    set(multiValueArgs DEPENDS WITH_REQUIREMENTS WITH)
    cmake_parse_arguments(
        PARSE_ARGV
        0
        ARGS
        "${options}"
        "${oneValueArgs}"
        "${multiValueArgs}"
    )

    if(NOT DEFINED ARGS_PYTHON_VERSION)
        set(ARGS_PYTHON_VERSION "3.13")
    endif()

    if(NOT DEFINED ARGS_PYTHON)
        message(SEND_ERROR "No python script specified")
        return()
    endif()

    if(NOT EXISTS ${ARGS_PYTHON})
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ARGS_PYTHON})
            set(ARGS_PYTHON ${CMAKE_CURRENT_SOURCE_DIR}/${ARGS_PYTHON})
        else()
            message(SEND_ERROR "Python script not found: ${ARGS_PYTHON}")
            return()
        endif()
    endif()

    set(_arg_isolated "")
    if(ARGS_ISOLATED)
        set(_arg_isolated "--isolated")
    endif()

    set(_depends "${ARGS_PYTHON}")
    set(_with_args "")
    foreach(_requirement ${ARGS_WITH_REQUIREMENTS})
        list(APPEND _depends ${_requirement})
        list(APPEND _with_args "--with-requirements;${_requirement}")
    endforeach()

    foreach(_requirement ${ARGS_WITH})
        list(APPEND _with_args "--with;${_requirement}")
        if(IS_DIRECTORY ${_requirement})
            if(NOT ACTS_USE_SYSTEM_LIBS)
                file(GLOB_RECURSE _depends_py ${_requirement}/*)
                list(APPEND _depends ${_depends_py})
            else()
                # If we are not using uv, then we use pip to install the
                # package into the virtual environment in the build directory.
                # The --no-build-isolation flag ensures that we don't
                # automatically download setuptools. The --no-index flag
                # ensures that nothing can be downloaded ever, and the
                # --no-deps is necessary to convince pip that the dependencies
                # are already in the environment.
                execute_process(
                    COMMAND
                        ${CMAKE_BINARY_DIR}/codegen_venv/bin/python -m pip
                        install --no-build-isolation --no-index --no-deps
                        ${_requirement}
                    OUTPUT_QUIET
                )
            endif()
        endif()
    endforeach()

    get_filename_component(_output_name ${ARGS_OUTPUT} NAME)

    string(SHA1 _output_hash ${_output_name})

    set(_codegen_root ${CMAKE_CURRENT_BINARY_DIR}/codegen/${_output_hash})
    set(_output_file ${_codegen_root}/${ARGS_OUTPUT})

    get_filename_component(_output_dir ${_output_file} DIRECTORY)
    file(MAKE_DIRECTORY ${_output_dir})

    if(NOT ACTS_USE_SYSTEM_LIBS)
        add_custom_command(
            OUTPUT ${_output_file}
            COMMAND
                env -i UV_NO_CACHE=1 ${uv_exe} run --quiet --python
                ${ARGS_PYTHON_VERSION} --no-project ${_arg_isolated}
                ${_with_args} ${ARGS_PYTHON} ${_output_file}
            DEPENDS ${_depends}
            COMMENT "Generating ${ARGS_OUTPUT}"
            VERBATIM
        )
    else()
        # If not using uv, just run Python from the virtual environment that
        # we created above.
        add_custom_command(
            OUTPUT ${_output_file}
            COMMAND
                ${CMAKE_BINARY_DIR}/codegen_venv/bin/python ${ARGS_PYTHON}
                ${_output_file}
            DEPENDS ${_depends}
            COMMENT "Generating ${ARGS_OUTPUT}"
            VERBATIM
        )
    endif()

    set(_internal_target codegen_${_output_hash}_Internal)
    add_custom_target(${_internal_target} DEPENDS ${_output_file})

    add_dependencies(${ARGS_ADD_TO_TARGET} ${_internal_target})
    target_include_directories(${ARGS_ADD_TO_TARGET} PRIVATE ${_codegen_root})

    # Add a central copde generation target that depends on all codegen targets, so that we can build only them in one go
    if(NOT TARGET ActsCodegen)
        add_custom_target(ActsCodegen)
    endif()
    add_dependencies(ActsCodegen ${_internal_target})
endfunction()
