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
  "include/traccc/cuda/utils/stream.hpp"
  "src/stream.cpp"
  "include/traccc/cuda/utils/make_magnetic_field.hpp"
  "src/make_magnetic_field.cpp"
  # Implementation details.
  "src/cuda_error_handling.hpp"
  "src/cuda_error_handling.cpp"
  "src/get_device.hpp"
  "src/get_device.cpp"
  "src/select_device.hpp"
  "src/select_device.cpp" )
target_link_libraries( traccc_cuda_utils
  PUBLIC traccc::core
  PRIVATE CUDA::cudart covfie::core covfie::cuda )
