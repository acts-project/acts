# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Find the HIP toolkit
find_package(HIPToolkit REQUIRED)

# Enable the HIP language, otherwise covfie::hip would cause an error.
enable_language(HIP)

# Project include(s).
include(traccc-compiler-options-cpp)

# Set up the build of the traccc::hip_utils library.
traccc_add_library( traccc_hip_utils hip_utils TYPE SHARED
  # Public facing code.
  "include/traccc/hip/utils/stream_wrapper.hpp"
  "src/stream_wrapper.cpp"
  "include/traccc/hip/utils/make_magnetic_field.hpp"
  "src/make_magnetic_field.cpp"
  # Implementation details.
  "src/hip_error_handling.hpp"
  "src/hip_error_handling.cpp"
)
target_link_libraries(
    traccc_hip_utils
    PUBLIC traccc::core
    PRIVATE HIP::hiprt covfie::core covfie::hip
)

# Check if hipStreamGetDevice(...) is available. It's done through a helper
# function, so that CMAKE_CXX_FLAGS "would be handled correctly".
include(CheckCXXSourceCompiles)

function(traccc_hip_utils_set_flags)
    # Needed not to pick up -Werror for the following test.
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_REQUIRED_LIBRARIES HIP::hiprt)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    check_cxx_source_compiles(
        "#include <hip/hip_runtime_api.h>
       void foo() {
          hipStream_t stream{nullptr};
          int device{-1};
          hipStreamGetDevice(stream, &device);
       }"
        TRACCC_HAVE_HIP_STREAM_GET_DEVICE
    )
    if(TRACCC_HAVE_HIP_STREAM_GET_DEVICE)
        target_compile_definitions(
            traccc_hip_utils
            PRIVATE TRACCC_HAVE_HIP_STREAM_GET_DEVICE
        )
    endif()
endfunction()

traccc_hip_utils_set_flags()
