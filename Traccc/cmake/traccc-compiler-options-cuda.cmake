# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# FindCUDAToolkit needs at least CMake 3.17, and C++17 support
# (set in the project's main CMakeLists.txt file) needs CMake 3.18.
cmake_minimum_required( VERSION 3.18 )

# Include the helper function(s).
include( traccc-functions )

# Figure out the properties of CUDA being used.
find_package( CUDAToolkit REQUIRED )

# Turn on the correct setting for the __cplusplus macro with MSVC.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   traccc_add_flag( CMAKE_CUDA_FLAGS "-Xcompiler /Zc:__cplusplus" )
endif()

if( "${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA" )
   traccc_add_flag( CMAKE_CUDA_FLAGS "-Wall" )
   traccc_add_flag( CMAKE_CUDA_FLAGS "-Wextra" )
   traccc_add_flag( CMAKE_CUDA_FLAGS "-Wconversion" )
endif()

# Allow to use functions in device code that are constexpr, even if they are
# not marked with __device__.
traccc_add_flag( CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr --use_fast_math" )

# Make CUDA generate debug symbols for the device code as well in a debug
# build.
traccc_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-G --keep" )

# Work around a bug in CUDA 12.8. Enabling the embedding of C++ source code in
# generated PTX code causes a ptxas error. A solution was promised for
# CUDA 13.1, but this has not yet surfaced.
#
# TODO: Add an upper bound to this statement when a fix in CUDA is presented.
if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "12.8")
   message(
      STATUS
      "Disabling C++ source in PTX in order to work around a bug in CUDA 12.8:"
   )
else()
   traccc_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-src-in-ptx" )
endif()

# Ensure that line information is embedded in debugging builds so that
# profilers have access to line data.
traccc_add_flag( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-lineinfo" )

# Fail on warnings, if asked for that behaviour.
if( TRACCC_FAIL_ON_WARNINGS )
   if( ( "${CUDAToolkit_VERSION}" VERSION_GREATER_EQUAL "10.2" ) AND
       ( "${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA" ) )
      traccc_add_flag( CMAKE_CUDA_FLAGS "-Werror all-warnings" )
   elseif( "${CMAKE_CUDA_COMPILER_ID}" MATCHES "Clang" )
      traccc_add_flag( CMAKE_CUDA_FLAGS "-Werror" )
   endif()
endif()
