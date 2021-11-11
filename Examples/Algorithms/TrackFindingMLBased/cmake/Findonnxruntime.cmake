# This will define the following variables:
#   onnxruntime_FOUND        -- True if the system has the onnxruntime library
#   onnxruntime_INCLUDE_DIRS -- The include directories for onnxruntime
#   onnxruntime_LIBRARIES    -- Libraries to link against

include(FindPackageHandleStandardArgs)

find_library(onnxruntime_LIBRARY
    NAMES onnxruntime
    HINTS ${onnxruntime_DIR}
    PATH_SUFFIXES lib lib32 lib64
    DOC "The onnxruntime library")

find_path(onnxruntime_INCLUDE_DIRS
    NAMES core/session/onnxruntime_cxx_api.h core/session/providers/cuda_provider_factory.h
    PATH_SUFFIXES include include/onnxruntime
    HITS ${onxxruntime_DIR}
    # PATHS /usr/local/include/onnxruntime/core/session
    # /usr/local/include/onnxruntime/core/session/providers 
    DOC "The onxxruntime include directory")

if (NOT onnxruntime_INCLUDE_DIRS)
    message(FATAL_ERROR "onnxruntime includes not found")
endif()

find_package_handle_standard_args(onnxruntime
    REQUIRED_VARS onnxruntime_LIBRARY onnxruntime_INCLUDE_DIRS)

# list(APPEND onnxruntime_INCLUDE_DIRS onnxruntime_API_DIR)
# list(APPEND onnxruntime_INCLUDE_DIRS onnxruntime_CUDA_DIR)

add_library(onnxruntime SHARED IMPORTED)
set_property(TARGET onnxruntime PROPERTY IMPORTED_LOCATION "${onnxruntime_LIBRARY}")
set_property(TARGET onnxruntime PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${onnxruntime_INCLUDE_DIRS}")

mark_as_advanced(onnxruntime_FOUND onnxruntime_INCLUDE_DIRS onnxruntime_LIBRARY)