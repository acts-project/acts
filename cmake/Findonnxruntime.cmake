# Find the ONNX Runtime include directory and library.
#
# This module defines the `onnxruntime` imported target that encodes all
# necessary information in its target properties.

find_library(
    onnxruntime_LIBRARY
    NAMES onnxruntime
    PATHS ${onnxruntime_DIR}
    PATH_SUFFIXES lib lib32 lib64
    DOC "The ONNXRuntime library"
)

if(NOT onnxruntime_LIBRARY)
    message(FATAL_ERROR "onnxruntime library not found")
else()
    message(STATUS "Found onnxruntime library at ${onnxruntime_LIBRARY}")
endif()

find_path(
    onnxruntime_INCLUDE_DIR
    NAMES onnxruntime_cxx_api.h
    PATHS ${onnxruntime_DIR}
    PATH_SUFFIXES
        include
        include/onnxruntime
        include/onnxruntime/core/session
        include/core/session
    DOC "The ONNXRuntime include directory"
)

if(NOT onnxruntime_INCLUDE_DIR)
    message(FATAL_ERROR "onnxruntime includes not found")
else()
    file(READ ${onnxruntime_INCLUDE_DIR}/onnxruntime_c_api.h ver)
    string(REGEX MATCH "ORT_API_VERSION ([0-9]*)" _ ${ver})
    set(onnxruntime_API_VERSION ${CMAKE_MATCH_1})
    message(
        STATUS
        "Found onnxruntime includes at ${onnxruntime_INCLUDE_DIR} (API version: ${onnxruntime_API_VERSION})"
    )
endif()

if("${onnxruntime_API_VERSION}" VERSION_LESS "${onnxruntime_MIN_API_VERSION}")
    message(
        FATAL_ERROR
        "onnxruntime API version ${onnxruntime_MIN_API_VERSION} or greater required"
    )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    onnxruntime
    REQUIRED_VARS onnxruntime_LIBRARY onnxruntime_INCLUDE_DIR
)

add_library(onnxruntime SHARED IMPORTED)
set_property(
    TARGET onnxruntime
    PROPERTY IMPORTED_LOCATION ${onnxruntime_LIBRARY}
)
set_property(
    TARGET onnxruntime
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${onnxruntime_INCLUDE_DIR}
)
set_property(
    TARGET onnxruntime
    PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${onnxruntime_INCLUDE_DIR}
)
mark_as_advanced(onnxruntime_FOUND onnxruntime_INCLUDE_DIR onnxruntime_LIBRARY)

add_library(onnxruntime::onnxruntime ALIAS onnxruntime)
