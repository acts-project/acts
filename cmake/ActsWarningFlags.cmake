# Compiler warning policy, shared by Acts, Detray and Traccc.
#
# This module is the single source of truth for *warning* flags and for
# warnings-as-errors in this repository. It is deliberately limited to that:
# functional flags (architectures, -march, --use_fast_math, -G/-lineinfo,
# --expt-relaxed-constexpr, ...) stay in the per-project compiler option files.
#
# Flags are applied with add_compile_options(), i.e. as a *directory* property,
# rather than by appending to the global CMAKE_<LANG>_FLAGS. That is what keeps
# them off the third-party code that Acts, Detray and Traccc pull in via
# add_subdirectory()/FetchContent: only directories that explicitly opt in by
# calling acts_apply_warning_flags() are affected.
#
# Detray and Traccc include this file by absolute path (via ../cmake) instead of
# through CMAKE_MODULE_PATH, so that Acts' Find<Package>.cmake modules do not
# shadow theirs.

include_guard(GLOBAL)

include(CheckCXXCompilerFlag)

# __COUNTER__ is classified as a C2y extension in Clang 22+. Third-party
# generated code (podio, Boost.Test, etc.) relies on it.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    check_cxx_compiler_flag("-Wno-c2y-extensions" ACTS_HAS_WNO_C2Y_EXTENSIONS)
endif()

# _acts_warning_flags(<out_var> <profile> <language>)
#
# Computes the warning flag list for one profile/language combination. This is
# the one place where the projects' warning sets are allowed to differ.
function(_acts_warning_flags out profile lang)
    set(flags "")

    if(profile STREQUAL "none")
        # nothing at all: explicit opt-out
    elseif(lang STREQUAL "CXX")
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
            # the common base, shared by all three projects
            list(
                APPEND flags
                -Wall
                -Wextra
                -Wpedantic
                -Wshadow
                -Wzero-as-null-pointer-constant
                -Wold-style-cast
                -Woverloaded-virtual
            )

            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                # useful conversion checks like float-to-bool, float-to-int.
                # Not enabled for GCC, where -Wfloat-conversion is much more
                # aggressive and also triggers on e.g. double-to-float.
                list(APPEND flags -Wfloat-conversion)
                if(ACTS_HAS_WNO_C2Y_EXTENSIONS)
                    list(APPEND flags -Wno-c2y-extensions)
                endif()
            endif()

            # -Wnull-dereference gets applied to -isystem includes in GCC13,
            # which causes lots of warnings in code we have no control over
            if(
                CMAKE_CXX_COMPILER_ID MATCHES "Clang"
                OR (
                    CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
                    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 12
                )
            )
                list(APPEND flags -Wnull-dereference)
            endif()

            # Detray and Traccc enforce a stricter set than Acts does. Acts'
            # own code is not -Wconversion clean, so it stays on the base set.
            if(profile STREQUAL "detray" OR profile STREQUAL "traccc")
                list(APPEND flags -Wconversion -Wunused-local-typedefs)
                if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                    list(APPEND flags -Wshorten-64-to-32)
                endif()
            endif()
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            list(APPEND flags /W4)
        endif()
    elseif(lang STREQUAL "CUDA")
        # Acts (Plugins/Gnn) and Detray do not warn on device code today.
        if(
            profile STREQUAL "traccc"
            AND CMAKE_CUDA_COMPILER_ID MATCHES "NVIDIA"
        )
            list(APPEND flags -Wall -Wextra -Wconversion)
        endif()
    elseif(lang STREQUAL "SYCL")
        if(
            CMAKE_SYCL_COMPILER_ID STREQUAL "IntelLLVM"
            OR CMAKE_SYCL_COMPILER_ID MATCHES "Clang"
        )
            list(
                APPEND flags
                -Wall
                -Wextra
                -Wno-unknown-cuda-version
                -Wshadow
                -Wunused-local-typedefs
            )
            if(NOT WIN32)
                list(APPEND flags -Wpedantic)
            endif()
            if(profile STREQUAL "traccc")
                list(APPEND flags -Wconversion)
            endif()
            if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
                list(APPEND flags -Wno-unused-command-line-argument)
            endif()
        endif()
    elseif(lang STREQUAL "HIP")
        if(CMAKE_HIP_PLATFORM MATCHES "^(amd|hcc)$")
            list(
                APPEND flags
                -Wall
                -Wextra
                -Wshadow
                -Wunused-local-typedefs
                -Wpedantic
            )
        endif()
    endif()

    set(${out} "${flags}" PARENT_SCOPE)
endfunction()

# acts_apply_warning_flags(
#     PROFILE           <acts|detray|traccc|none>
#     LANGUAGES         <CXX|CUDA|SYCL|HIP>...
#     [WARNING_AS_ERROR <bool>]   # defaults to the profile's project option
# )
#
# Applies the warning flags for PROFILE to the current directory and everything
# added below it, and enables warnings-as-errors in the same scope if the
# project's <PROJECT>_FAIL_ON_WARNINGS option is set.
#
# This must be called before any target is created in the directory and before
# any add_subdirectory(): both the directory COMPILE_OPTIONS property and the
# CMAKE_COMPILE_WARNING_AS_ERROR variable are snapshotted at those points.
#
# This is a macro rather than a function because CMAKE_COMPILE_WARNING_AS_ERROR
# and the MSVC /W<n> fixup have to land in the caller's directory scope.
macro(acts_apply_warning_flags)
    cmake_parse_arguments(
        _awf
        ""
        "PROFILE;WARNING_AS_ERROR"
        "LANGUAGES"
        ${ARGN}
    )

    if(NOT _awf_PROFILE)
        message(FATAL_ERROR "acts_apply_warning_flags: PROFILE is required")
    elseif(NOT _awf_PROFILE MATCHES "^(acts|detray|traccc|none)$")
        message(
            FATAL_ERROR
            "acts_apply_warning_flags: unknown PROFILE '${_awf_PROFILE}'"
        )
    endif()
    if(NOT _awf_LANGUAGES)
        message(FATAL_ERROR "acts_apply_warning_flags: LANGUAGES is required")
    endif()

    # per-project switches
    set(_awf_enabled ON)
    set(_awf_wae OFF)
    if(_awf_PROFILE STREQUAL "acts")
        set(_awf_enabled ${ACTS_ENABLE_WARNINGS})
        set(_awf_wae ${ACTS_FAIL_ON_WARNINGS})
    elseif(_awf_PROFILE STREQUAL "detray")
        set(_awf_wae ${DETRAY_FAIL_ON_WARNINGS})
    elseif(_awf_PROFILE STREQUAL "traccc")
        set(_awf_wae ${TRACCC_FAIL_ON_WARNINGS})
    elseif(_awf_PROFILE STREQUAL "none")
        set(_awf_enabled OFF)
    endif()
    if(DEFINED _awf_WARNING_AS_ERROR)
        set(_awf_wae ${_awf_WARNING_AS_ERROR})
    endif()
    if(NOT _awf_enabled)
        set(_awf_PROFILE "none")
        set(_awf_wae OFF)
    endif()

    foreach(_awf_lang IN LISTS _awf_LANGUAGES)
        _acts_warning_flags(_awf_flags "${_awf_PROFILE}" "${_awf_lang}")
        if(_awf_flags)
            if(
                _awf_lang STREQUAL "CXX"
                AND CMAKE_CXX_COMPILER_ID MATCHES "MSVC"
            )
                # avoid D9025 ("overriding /W3 with /W4") on every TU
                string(
                    REGEX REPLACE "/W[0-9]"
                    ""
                    CMAKE_CXX_FLAGS
                    "${CMAKE_CXX_FLAGS}"
                )
            endif()
            add_compile_options(
                "$<$<COMPILE_LANGUAGE:${_awf_lang}>:${_awf_flags}>"
            )
        endif()

        # announce each profile once, not once per call site
        get_property(
            _awf_announced
            GLOBAL
            PROPERTY "_acts_warnings_${_awf_PROFILE}_${_awf_lang}"
        )
        if(NOT _awf_announced)
            message(
                STATUS
                "Warning flags [${_awf_PROFILE}/${_awf_lang}]: ${_awf_flags}"
            )
            set_property(
                GLOBAL
                PROPERTY "_acts_warnings_${_awf_PROFILE}_${_awf_lang}" ON
            )
        endif()
    endforeach()

    if(_awf_wae)
        if(CMAKE_VERSION VERSION_LESS 3.24)
            message(
                FATAL_ERROR
                "Warnings-as-errors requires CMake >= 3.24, but this is "
                "CMake ${CMAKE_VERSION}. Configure without "
                "<PROJECT>_FAIL_ON_WARNINGS, or use a newer CMake."
            )
        endif()
        # The native mechanism picks the right flag per compiler (-Werror,
        # /WX, or `-Werror all-warnings` for nvcc) and keeps
        # `cmake --build ... --compile-no-warning-as-error` working.
        set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

        # ... except for SYCL, which is a vecmem-provided pseudo-language that
        # CMake has no COMPILE_WARNING_AS_ERROR implementation for. The
        # property would be silently ignored there.
        if("SYCL" IN_LIST _awf_LANGUAGES AND NOT _awf_PROFILE STREQUAL "none")
            add_compile_options($<$<COMPILE_LANGUAGE:SYCL>:-Werror>)
        endif()
    endif()

    unset(_awf_PROFILE)
    unset(_awf_LANGUAGES)
    unset(_awf_WARNING_AS_ERROR)
    unset(_awf_UNPARSED_ARGUMENTS)
    unset(_awf_KEYWORDS_MISSING_VALUES)
    unset(_awf_enabled)
    unset(_awf_wae)
    unset(_awf_lang)
    unset(_awf_flags)
    unset(_awf_announced)
endmacro()
