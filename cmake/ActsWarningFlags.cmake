# Shared warning flags for Acts, Detray and Traccc.
# Apply per directory after dependency setup to keep warnings off third-party code.
# Other compiler options stay in each project's compiler-options files.

include_guard(GLOBAL)

include(CheckCXXCompilerFlag)

# Allow __COUNTER__ in third-party headers with Clang 22+.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    check_cxx_compiler_flag("-Wno-c2y-extensions" ACTS_HAS_WNO_C2Y_EXTENSIONS)
endif()

# _acts_warning_flags(<out_var> <profile> <language>)

function(_acts_warning_flags out profile lang)
    set(flags "")

    if(profile STREQUAL "none")
        # No warning flags.
    elseif(lang STREQUAL "CXX")
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
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
                # GCC also warns on double-to-float conversions, which Acts allows.
                list(APPEND flags -Wfloat-conversion)
                if(ACTS_HAS_WNO_C2Y_EXTENSIONS)
                    list(APPEND flags -Wno-c2y-extensions)
                endif()
            endif()

            # GCC 13+ reports null dereferences in system headers.
            if(
                CMAKE_CXX_COMPILER_ID MATCHES "Clang"
                OR (
                    CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
                    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 12
                )
            )
                list(APPEND flags -Wnull-dereference)
            endif()

            # Acts is not yet -Wconversion clean.
            if(profile STREQUAL "detray" OR profile STREQUAL "traccc")
                list(APPEND flags -Wconversion -Wunused-local-typedefs)
                # Suppress GCC 13+ false positives in Eigen for Detray.
                if(
                    profile STREQUAL "detray"
                    AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
                    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 12
                )
                    list(APPEND flags -Wno-maybe-uninitialized)
                endif()
                if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                    list(APPEND flags -Wshorten-64-to-32)
                endif()
            endif()
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            list(APPEND flags /W4)
        endif()
    elseif(lang STREQUAL "CUDA")
        # Only Traccc enables additional CUDA warnings.
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
# Call after dependency setup and before adding project targets/subdirectories.
# Uses <PROJECT>_FAIL_ON_WARNINGS to control warnings-as-errors.
# A macro is needed to set compiler variables in the caller's directory scope.
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

    set(_awf_wae OFF)
    if(_awf_PROFILE STREQUAL "acts")
        set(_awf_wae ${ACTS_FAIL_ON_WARNINGS})
    elseif(_awf_PROFILE STREQUAL "detray")
        set(_awf_wae ${DETRAY_FAIL_ON_WARNINGS})
    elseif(_awf_PROFILE STREQUAL "traccc")
        set(_awf_wae ${TRACCC_FAIL_ON_WARNINGS})
    endif()
    if(DEFINED _awf_WARNING_AS_ERROR)
        set(_awf_wae ${_awf_WARNING_AS_ERROR})
    endif()
    if(_awf_PROFILE STREQUAL "none")
        set(_awf_wae OFF)
    endif()

    foreach(_awf_lang IN LISTS _awf_LANGUAGES)
        _acts_warning_flags(_awf_flags "${_awf_PROFILE}" "${_awf_lang}")
        if(_awf_flags)
            if(
                _awf_lang STREQUAL "CXX"
                AND CMAKE_CXX_COMPILER_ID MATCHES "MSVC"
            )
                # Avoid MSVC D9025 when replacing /W3 with /W4.
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

        # Log each profile/language pair once.
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
        # Let CMake select the compiler flag and honor --compile-no-warning-as-error.
        set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

        # CMake does not implement COMPILE_WARNING_AS_ERROR for vecmem's SYCL language.
        if("SYCL" IN_LIST _awf_LANGUAGES AND NOT _awf_PROFILE STREQUAL "none")
            add_compile_options($<$<COMPILE_LANGUAGE:SYCL>:-Werror>)
        endif()
    endif()

    unset(_awf_PROFILE)
    unset(_awf_LANGUAGES)
    unset(_awf_WARNING_AS_ERROR)
    unset(_awf_UNPARSED_ARGUMENTS)
    unset(_awf_KEYWORDS_MISSING_VALUES)
    unset(_awf_wae)
    unset(_awf_lang)
    unset(_awf_flags)
    unset(_awf_announced)
endmacro()
