# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# FindCUDAToolkit needs at least CMake 3.17, and C++17 support
# (set in the project's main CMakeLists.txt file) needs CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Include the helper function(s).
include(traccc-functions)

# Figure out the properties of CUDA being used.
find_package(CUDAToolkit REQUIRED)

# Turn on the correct setting for the __cplusplus macro with MSVC.
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    traccc_add_flag( CMAKE_CUDA_FLAGS "-Xcompiler /Zc:__cplusplus" )
endif()

# Allow to use functions in device code that are constexpr, even if they are
# not marked with __device__.
traccc_add_flag( CMAKE_CUDA_FLAGS "--use_fast_math" )

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

# Warnings and warnings-as-errors come from the shared ActsWarningFlags.cmake
# module, and are applied to the directory that includes this file. CMake picks
# the right flag for the compiler in use (`-Werror all-warnings` for
# nvcc >= 10.2, `-Werror` for clang-cuda).
acts_apply_warning_flags(PROFILE traccc LANGUAGES CUDA)
