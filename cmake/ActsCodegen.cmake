include_guard(GLOBAL)

# Ensure this module's own directory is on the module path so that nested bare
# includes (e.g. ActsEnsureUv) resolve even when this file is included by
# absolute path from outside the ACTS module path — for example detray's build,
# which manages its own CMAKE_MODULE_PATH. include_guard already runs this once;
# the IN_LIST check just avoids a redundant entry in the normal ACTS build.
if(NOT "${CMAKE_CURRENT_LIST_DIR}" IN_LIST CMAKE_MODULE_PATH)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
endif()

# Root of the ACTS source tree; this file lives in <root>/cmake. This is the
# anchor for the pre-generated code layout below: it resolves to the same
# directory for the ACTS build and for detray, whether detray is built via
# add_subdirectory from ACTS or standalone inside the ACTS tree. Both therefore
# agree on where a given generated file lives inside a source archive.
#
# Everything this module sets outside of a function has to live in the cache:
# the include guard is global, so the plain variables would only ever exist in
# the scope of whichever directory happened to include this file first.
get_filename_component(
    _acts_codegen_source_root
    "${CMAKE_CURRENT_LIST_DIR}/.."
    ABSOLUTE
)
set(ACTS_CODEGEN_SOURCE_ROOT "${_acts_codegen_source_root}" CACHE INTERNAL "")

# The manifest is the single source of truth for what gets generated: this
# module drives acts_code_generation from it, and CI/pregenerate_codegen.py
# produces files from the same list, so neither side can invent a unit or a
# location the other does not know about.
file(
    READ "${ACTS_CODEGEN_SOURCE_ROOT}/codegen/manifest.json"
    _acts_codegen_manifest
)

# A global property rather than a cache entry: cache values must not contain
# newlines, and CMake silently truncates them if they do.
set_property(GLOBAL PROPERTY ACTS_CODEGEN_MANIFEST "${_acts_codegen_manifest}")

# Accessor, so that directories which only see this module through the global
# include guard -- and the CI fixture in CI/codegen_prebuilt_check -- can read
# the manifest without knowing where it lives.
function(acts_codegen_manifest out_var)
    get_property(_manifest GLOBAL PROPERTY ACTS_CODEGEN_MANIFEST)
    set(${out_var} "${_manifest}" PARENT_SCOPE)
endfunction()

# Read a JSON array of source-root-relative paths out of a manifest unit and
# return it as a list of absolute paths.
function(_acts_codegen_manifest_paths out_var unit member)
    set(_result "")
    string(JSON _count LENGTH "${unit}" ${member})
    if(_count GREATER 0)
        math(EXPR _last "${_count} - 1")
        foreach(_i RANGE ${_last})
            string(JSON _entry GET "${unit}" ${member} ${_i})
            list(APPEND _result "${ACTS_CODEGEN_SOURCE_ROOT}/${_entry}")
        endforeach()
    endif()
    set(${out_var} "${_result}" PARENT_SCOPE)
endfunction()

# Source archives built by CI/make_source_archive.sh ship the generated code in
# this directory below the source root. When it is present, the generators are
# not run at all — that is what lets a build from such an archive work without
# network access, uv, or a Python environment.
set(ACTS_CODEGEN_PREBUILT_DIRNAME "prebuilt-codegen" CACHE INTERNAL "")

# ACTS_CODEGEN_PREBUILT_DIR and ACTS_CODEGEN_REQUIRE_PREBUILT are declared in the
# top-level CMakeLists.txt, where all options live. They are only read here, so
# that this module keeps working when detray includes it standalone and neither
# of them is defined.
#
# Resolve the directory to actually use: what the option says if it was set, and
# otherwise the one a source archive ships.
if(ACTS_CODEGEN_PREBUILT_DIR)
    set(_acts_codegen_prebuilt "${ACTS_CODEGEN_PREBUILT_DIR}")
elseif(EXISTS "${ACTS_CODEGEN_SOURCE_ROOT}/${ACTS_CODEGEN_PREBUILT_DIRNAME}")
    set(_acts_codegen_prebuilt
        "${ACTS_CODEGEN_SOURCE_ROOT}/${ACTS_CODEGEN_PREBUILT_DIRNAME}"
    )
else()
    set(_acts_codegen_prebuilt "")
endif()
set(ACTS_CODEGEN_PREBUILT_RESOLVED
    "${_acts_codegen_prebuilt}"
    CACHE INTERNAL
    ""
)

if(ACTS_CODEGEN_PREBUILT_RESOLVED)
    message(
        STATUS
        "Codegen: using pre-generated code from ${ACTS_CODEGEN_PREBUILT_RESOLVED}"
    )
elseif(ACTS_CODEGEN_REQUIRE_PREBUILT)
    # Most likely cause: built from the source archive GitHub generates for a
    # tag, which does not carry the generated code, rather than from the archive
    # attached to the release.
    message(
        FATAL_ERROR
        "ACTS_CODEGEN_REQUIRE_PREBUILT is set, but no pre-generated code was found in "
        "${ACTS_CODEGEN_SOURCE_ROOT}/${ACTS_CODEGEN_PREBUILT_DIRNAME}. Source archives "
        "attached to an ACTS release contain that directory; the archives GitHub "
        "generates automatically for a tag do not. Either build from the release "
        "asset, point ACTS_CODEGEN_PREBUILT_DIR at a directory of pre-generated code, "
        "or turn ACTS_CODEGEN_REQUIRE_PREBUILT off to run the generators."
    )
endif()

# Set up whatever is needed to actually *run* the generators: either uv (which
# downloads itself, a Python interpreter and the generator dependencies) or a
# virtual environment derived from the ambient Python.
#
# This is deliberately lazy and only called from acts_code_generation when a
# generated file has no pre-generated counterpart. A build that is fully covered
# by the pre-generated code never calls it, and so never downloads anything
# and never needs a Python environment.
function(_acts_codegen_ensure_runtime)
    get_property(_ready GLOBAL PROPERTY ACTS_CODEGEN_RUNTIME_READY)
    if(_ready)
        return()
    endif()
    set_property(GLOBAL PROPERTY ACTS_CODEGEN_RUNTIME_READY TRUE)

    if(NOT ACTS_USE_SYSTEM_LIBS)
        message(STATUS "Configuring codegen")

        if(NOT DEFINED ACTS_CODEGEN_TMPDIR OR ACTS_CODEGEN_TMPDIR STREQUAL "")
            find_program(MKTEMP_EXE NAMES mktemp REQUIRED)
            execute_process(
                COMMAND ${MKTEMP_EXE} -d
                OUTPUT_VARIABLE _acts_codegen_tmpdir
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            set(ACTS_CODEGEN_TMPDIR
                "${_acts_codegen_tmpdir}"
                CACHE PATH
                "Codegen temporary directory for ACTS code generation"
            )
        endif()

        file(MAKE_DIRECTORY "${ACTS_CODEGEN_TMPDIR}")
        message(STATUS "Codegen temporary directory: ${ACTS_CODEGEN_TMPDIR}")

        include(ActsEnsureUv)
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
            COMMAND
                ${Python_EXECUTABLE} -m venv ${CMAKE_BINARY_DIR}/codegen_venv
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
            WRITE "${_python_nested_package_dir}/_base_packages.pth"
            ${_python_package_dir}
        )

        message(
            STATUS
            "Virtual environment based on ${_python_package_dir} created in ${CMAKE_BINARY_DIR}/codegen_venv/"
        )
    endif()
endfunction()

function(acts_code_generation)
    set(oneValueArgs ADD_TO_TARGET KEY RESULT_INCLUDE_DIR RESULT_OUTPUT)
    cmake_parse_arguments(PARSE_ARGV 0 ARGS "" "${oneValueArgs}" "")

    if(NOT DEFINED ARGS_KEY)
        message(SEND_ERROR "acts_code_generation: no KEY specified")
        return()
    endif()
    set(_key "${ARGS_KEY}")

    # Everything about a generated file other than which target consumes it
    # comes from the manifest, so that this build and whoever pre-generates the
    # code cannot end up with different ideas about what is generated or where
    # it goes.
    acts_codegen_manifest(_manifest)
    string(
        JSON _unit
        ERROR_VARIABLE _unit_error
        GET "${_manifest}"
        units
        "${_key}"
    )
    if(_unit_error)
        message(
            SEND_ERROR
            "acts_code_generation: '${_key}' is not in ${ACTS_CODEGEN_SOURCE_ROOT}/codegen/manifest.json. "
            "Add it there rather than passing the details here, otherwise pre-generated code will not cover it."
        )
        return()
    endif()

    string(JSON _script GET "${_unit}" script)
    string(JSON ARGS_OUTPUT GET "${_unit}" output)
    string(JSON ARGS_PYTHON_VERSION GET "${_unit}" python_version)
    string(JSON _isolated GET "${_unit}" isolated)

    set(ARGS_PYTHON "${ACTS_CODEGEN_SOURCE_ROOT}/${_script}")
    if(NOT EXISTS ${ARGS_PYTHON})
        message(SEND_ERROR "Python script not found: ${ARGS_PYTHON}")
        return()
    endif()

    _acts_codegen_manifest_paths(
        ARGS_WITH_REQUIREMENTS
        "${_unit}"
        with_requirements
    )
    _acts_codegen_manifest_paths(ARGS_WITH "${_unit}" with)

    get_filename_component(_output_name ${ARGS_OUTPUT} NAME)

    string(SHA1 _output_hash ${_output_name})

    set(_codegen_root ${CMAKE_CURRENT_BINARY_DIR}/codegen/${_output_hash})
    set(_output_file ${_codegen_root}/${ARGS_OUTPUT})

    get_filename_component(_output_dir ${_output_file} DIRECTORY)
    file(MAKE_DIRECTORY ${_output_dir})

    set(_prebuilt "")
    if(
        ACTS_CODEGEN_PREBUILT_RESOLVED
        AND EXISTS "${ACTS_CODEGEN_PREBUILT_RESOLVED}/${_key}"
    )
        set(_prebuilt "${ACTS_CODEGEN_PREBUILT_RESOLVED}/${_key}")
    endif()

    if(_prebuilt)
        # Nothing to generate: copy the shipped file into the place the build
        # expects it. Everything below this branch is unchanged, so the target
        # wiring does not care where the file came from.
        add_custom_command(
            OUTPUT ${_output_file}
            COMMAND
                ${CMAKE_COMMAND} -E copy_if_different ${_prebuilt}
                ${_output_file}
            DEPENDS ${_prebuilt}
            COMMENT "Using pre-generated ${ARGS_OUTPUT}"
            VERBATIM
        )
    else()
        if(ACTS_CODEGEN_REQUIRE_PREBUILT)
            message(
                FATAL_ERROR
                "No pre-generated file for ${_key} in ${ACTS_CODEGEN_PREBUILT_RESOLVED}. "
                "The pre-generated code does not cover this configuration, and "
                "ACTS_CODEGEN_REQUIRE_PREBUILT forbids running the generator."
            )
        elseif(ACTS_CODEGEN_PREBUILT_RESOLVED)
            message(
                WARNING
                "No pre-generated file for ${_key} in ${ACTS_CODEGEN_PREBUILT_RESOLVED}; "
                "falling back to running the generator, which needs uv or a Python environment"
            )
        endif()

        _acts_codegen_ensure_runtime()

        set(_arg_isolated "")
        if(_isolated)
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

        if(NOT ACTS_USE_SYSTEM_LIBS)
            # If using uv, run it in a clean environment but include the variables
            # that specify a proxy.
            set(_propagate
                HTTP_PROXY
                HTTPS_PROXY
                ALL_PROXY
                NO_PROXY
                SSL_CERT_FILE
            )
            foreach(_var IN LISTS _propagate)
                if(DEFINED ENV{${_var}})
                    list(APPEND _uv_environment "${_var}=$ENV{${_var}}")
                endif()
            endforeach()
            add_custom_command(
                OUTPUT ${_output_file}
                COMMAND
                    env -i UV_NO_CACHE=1
                    UV_PYTHON_INSTALL_DIR=${ACTS_CODEGEN_TMPDIR}/python_install_dir
                    ${_uv_environment} ${uv_exe} run --quiet --python
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
    endif()

    set(_internal_target codegen_${_output_hash}_Internal)
    add_custom_target(${_internal_target} DEPENDS ${_output_file})

    add_dependencies(${ARGS_ADD_TO_TARGET} ${_internal_target})

    # INTERFACE libraries (e.g. header-only detray::core) cannot carry PRIVATE
    # include directories; expose the generated headers on their INTERFACE so
    # consumers can find them. Compiled targets keep using PRIVATE as before.
    get_target_property(_add_to_target_type ${ARGS_ADD_TO_TARGET} TYPE)
    if(_add_to_target_type STREQUAL "INTERFACE_LIBRARY")
        target_include_directories(
            ${ARGS_ADD_TO_TARGET}
            INTERFACE $<BUILD_INTERFACE:${_codegen_root}>
        )

        # Use the CMake 3.19 functionality of attaching private sources to
        # interface targets to create an explicit dependency between the target
        # and the output file.
        target_sources(${ARGS_ADD_TO_TARGET} PRIVATE ${_output_file})
    else()
        target_include_directories(
            ${ARGS_ADD_TO_TARGET}
            PRIVATE ${_codegen_root}
        )
    endif()

    # Optionally expose the generated include root to the caller (used to set up
    # installation of the generated headers).
    if(DEFINED ARGS_RESULT_INCLUDE_DIR)
        set(${ARGS_RESULT_INCLUDE_DIR} ${_codegen_root} PARENT_SCOPE)
    endif()

    # The output path is only known here now that it comes from the manifest;
    # hand it back for callers that need to install the generated file.
    if(DEFINED ARGS_RESULT_OUTPUT)
        set(${ARGS_RESULT_OUTPUT} ${ARGS_OUTPUT} PARENT_SCOPE)
    endif()

    # Add a central copde generation target that depends on all codegen targets, so that we can build only them in one go
    if(NOT TARGET ActsCodegen)
        add_custom_target(ActsCodegen)
    endif()
    add_dependencies(ActsCodegen ${_internal_target})
endfunction()
