# Shared warning flags for Acts, Detray and Traccc.
# Apply per directory after dependency setup to keep warnings off third-party code.
# Other compiler options stay in each project's compiler-options files.

include_guard(GLOBAL)

include(CheckCXXCompilerFlag)

# Allow __COUNTER__ in third-party headers with Clang 22+.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    check_cxx_compiler_flag("-Wno-c2y-extensions" ACTS_HAS_WNO_C2Y_EXTENSIONS)
endif()

# Apply the same warnings to every project, scoped to the source language.
# Call after dependency setup, before creating the project's targets.
# A macro keeps CMAKE_COMPILE_WARNING_AS_ERROR in the caller's scope.
macro(acts_apply_warning_flags lang warnings_as_errors)
    set(_acts_warning_flags "")
    if("${lang}" STREQUAL "CXX")
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
            set(_acts_warning_flags
                -Wall
                -Wextra
                -Wpedantic
                -Wshadow
                -Wzero-as-null-pointer-constant
                -Wold-style-cast
                -Woverloaded-virtual
                -Wunused-local-typedefs
            )
            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                list(APPEND _acts_warning_flags -Wfloat-conversion)
                if(ACTS_HAS_WNO_C2Y_EXTENSIONS)
                    list(APPEND _acts_warning_flags -Wno-c2y-extensions)
                endif()
            endif()
            # GCC 13+ reports null dereferences in system headers.
            if(
                CMAKE_CXX_COMPILER_ID MATCHES "Clang"
                OR (
                    CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
                    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13
                )
            )
                list(APPEND _acts_warning_flags -Wnull-dereference)
            endif()
        elseif(MSVC)
            string(
                REGEX REPLACE
                "/W[0-9]"
                ""
                CMAKE_CXX_FLAGS
                "${CMAKE_CXX_FLAGS}"
            )
            list(APPEND _acts_warning_flags /W4)
        endif()
    elseif("${lang}" STREQUAL "CUDA")
        if(CMAKE_CUDA_COMPILER_ID MATCHES "NVIDIA")
            list(APPEND _acts_warning_flags -Wall -Wextra)
        endif()
    elseif("${lang}" STREQUAL "SYCL")
        if(
            CMAKE_SYCL_COMPILER_ID STREQUAL "IntelLLVM"
            OR CMAKE_SYCL_COMPILER_ID MATCHES "Clang"
        )
            list(
                APPEND
                _acts_warning_flags
                -Wall
                -Wextra
                -Wno-unknown-cuda-version
                -Wshadow
                -Wunused-local-typedefs
            )
            if(NOT WIN32)
                list(APPEND _acts_warning_flags -Wpedantic)
            endif()
            if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
                list(
                    APPEND
                    _acts_warning_flags
                    -Wno-unused-command-line-argument
                )
            endif()
        endif()
    elseif("${lang}" STREQUAL "HIP")
        if(CMAKE_HIP_PLATFORM MATCHES "^(amd|hcc)$")
            list(
                APPEND
                _acts_warning_flags
                -Wall
                -Wextra
                -Wshadow
                -Wunused-local-typedefs
                -Wpedantic
            )
        endif()
    endif()

    add_compile_options("$<$<COMPILE_LANGUAGE:${lang}>:${_acts_warning_flags}>")
    if(${warnings_as_errors})
        if(CMAKE_VERSION VERSION_LESS 3.24)
            message(FATAL_ERROR "Warnings-as-errors requires CMake >= 3.24")
        endif()
        # CMake selects the compiler flag and honors --compile-no-warning-as-error.
        set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
        # vecmem's custom SYCL language does not implement this CMake property.
        if("${lang}" STREQUAL "SYCL")
            add_compile_options("$<$<COMPILE_LANGUAGE:SYCL>:-Werror>")
        endif()
    endif()
    unset(_acts_warning_flags)
endmacro()
