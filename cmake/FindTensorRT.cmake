# ~~~
# Locate the TensorRT runtime.
#
# ACTS only ever deserializes pre-built ``.engine`` files, so this module
# deliberately looks for the runtime libraries only. The graph parsers
# (``nvonnxparser``, ``nvparsers``) are not searched for: building an engine
# from ONNX is an offline step performed outside of ACTS, and ``nvparsers`` no
# longer ships as of TensorRT 10.
#
# This module defines the following variables:
#
# - TensorRT_FOUND: A boolean specifying whether or not TensorRT was found.
# - TensorRT_VERSION: The exact version of TensorRT found.
# - TensorRT_INCLUDE_DIRS: The path to the TensorRT ``include`` folder.
# - TensorRT_LIBRARY_DIRS: The path to the TensorRT library directory.
#
# This module creates the following targets:
#
# - trt::nvinfer
# - trt::nvinfer_plugin
#
# Hints
# ^^^^^
# A user may set ``TensorRT_ROOT`` to an installation root to tell this module
# where to look.
# ~~~

find_path(
    TensorRT_INCLUDE_DIR
    NAMES NvInfer.h
    HINTS ${TensorRT_ROOT}
    ENV TensorRT_ROOT
    PATH_SUFFIXES include
)

if(EXISTS "${TensorRT_INCLUDE_DIR}/NvInferVersion.h")
    file(READ "${TensorRT_INCLUDE_DIR}/NvInferVersion.h" _trt_version_header)
    set(_trt_version_parts)
    foreach(_part MAJOR MINOR PATCH BUILD)
        string(
            REGEX MATCH "NV_TENSORRT_${_part}[ \t]+([0-9]+)"
            _trt_version_match
            "${_trt_version_header}"
        )
        list(APPEND _trt_version_parts "${CMAKE_MATCH_1}")
    endforeach()
    list(JOIN _trt_version_parts "." TensorRT_VERSION)
    unset(_trt_version_header)
    unset(_trt_version_parts)
    unset(_trt_version_match)
endif()

foreach(_component nvinfer nvinfer_plugin)
    find_library(
        TensorRT_${_component}_LIBRARY
        NAMES ${_component}
        HINTS ${TensorRT_ROOT}
        ENV TensorRT_ROOT
        PATH_SUFFIXES lib lib64
    )
    mark_as_advanced(TensorRT_${_component}_LIBRARY)
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    TensorRT
    REQUIRED_VARS
        TensorRT_INCLUDE_DIR
        TensorRT_nvinfer_LIBRARY
        TensorRT_nvinfer_plugin_LIBRARY
    VERSION_VAR TensorRT_VERSION
)

if(TensorRT_FOUND)
    set(TensorRT_INCLUDE_DIRS "${TensorRT_INCLUDE_DIR}")

    # Cached because ActsCreateSetup substitutes this into the setup scripts from
    # the top-level scope, while find_package(TensorRT) runs down in Plugins/Gnn.
    get_filename_component(
        _trt_library_dir
        "${TensorRT_nvinfer_LIBRARY}"
        DIRECTORY
    )
    set(TensorRT_LIBRARY_DIRS
        "${_trt_library_dir}"
        CACHE INTERNAL
        "TensorRT library directory"
    )
    unset(_trt_library_dir)

    foreach(_component nvinfer nvinfer_plugin)
        if(NOT TARGET trt::${_component})
            # Imported targets treat their interface includes as SYSTEM by
            # default, which keeps the TensorRT headers from tripping -Werror.
            add_library(trt::${_component} UNKNOWN IMPORTED)
            set_target_properties(
                trt::${_component}
                PROPERTIES
                    IMPORTED_LOCATION "${TensorRT_${_component}_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${TensorRT_INCLUDE_DIRS}"
            )
        endif()
    endforeach()
endif()

mark_as_advanced(TensorRT_INCLUDE_DIR)
