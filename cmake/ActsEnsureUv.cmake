include_guard(GLOBAL)

find_program(uv_exe uv)

if(uv_exe STREQUAL "uv_exe-NOTFOUND")
    message(STATUS "uv not found, installing it")

    set(_uv_version "0.7.19")
    set(_base_url
        "https://github.com/astral-sh/uv/releases/download/${_uv_version}"
    )

    if(NOT APPLE AND NOT UNIX)
        message(FATAL_ERROR "Unsupported platform: ${CMAKE_SYSTEM_NAME}")
    endif()

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
        if(APPLE)
            set(UV_URL "${_base_url}/uv-x86_64-apple-darwin.tar.gz")
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

execute_process(
    COMMAND ${uv_exe} --version
    OUTPUT_VARIABLE _uv_version_str
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Found uv ${_uv_version_str}: ${uv_exe}")
