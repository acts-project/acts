# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Find the CUDA toolkit
find_package( CUDAToolkit REQUIRED )

# Enable the CUDA language, otherwise covfie::cuda would cause an error.
enable_language( CUDA )

# Project include(s).
include( traccc-compiler-options-cpp )

# Set up the build of the traccc::cuda_utils library.
traccc_add_library( traccc_cuda_utils cuda_utils TYPE SHARED
  # Public facing code.
  "include/traccc/cuda/utils/stream_wrapper.hpp"
  "src/stream_wrapper.cpp"
  "include/traccc/cuda/utils/make_magnetic_field.hpp"
  "src/make_magnetic_field.cpp"
  # Implementation details.
  "src/cuda_error_handling.hpp"
  "src/cuda_error_handling.cpp" )
target_link_libraries( traccc_cuda_utils
  PUBLIC traccc::core
  PRIVATE CUDA::cudart covfie::core covfie::cuda )

# Check if cudaStreamGetDevice(...) is available. It's done through a helper
# function, so that CMAKE_CXX_FLAGS "would be handled correctly".
include(CheckCXXSourceCompiles)

function(traccc_cuda_utils_set_flags)
  set(CMAKE_CXX_FLAGS "") # Needed not to pick up -Werror for the following
                          # test.
  set(CMAKE_REQUIRED_INCLUDES "${CUDAToolkit_INCLUDE_DIRS}")
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  check_cxx_source_compiles(
    "#include <cuda_runtime_api.h>
     void foo() {
       cudaStream_t stream{nullptr};
       int device{-1};
       cudaStreamGetDevice(stream, &device);
     }"
    TRACCC_HAVE_CUDA_STREAM_GET_DEVICE)
  if(TRACCC_HAVE_CUDA_STREAM_GET_DEVICE)
    target_compile_definitions(traccc_cuda_utils
      PRIVATE TRACCC_HAVE_CUDA_STREAM_GET_DEVICE)
  endif()
endfunction()

traccc_cuda_utils_set_flags()
