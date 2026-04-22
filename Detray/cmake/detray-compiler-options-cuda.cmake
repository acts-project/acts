# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.21)

# Only set these compiler flags if we are the top level project.
if(PROJECT_IS_TOP_LEVEL)
    # Include the helper function(s).
    include(detray-functions)

    # Figure out the properties of CUDA being used.
    find_package(CUDAToolkit REQUIRED)

    # Turn on the correct setting for the __cplusplus macro with MSVC.
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
        detray_add_flag( CMAKE_CUDA_FLAGS "-Xcompiler /Zc:__cplusplus" )
    endif()

    # Set the CUDA architecture to build code for.
    set(CMAKE_CUDA_ARCHITECTURES
        "52"
        CACHE STRING
        "CUDA architectures to build device code for"
    )

    if("${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA")
        # Allow to use functions in device code that are constexpr, even if they are
        # not marked with __device__.
        detray_add_flag( CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr" )
    endif()

    # Make CUDA generate debug symbols for the device code as well in a debug
    # build.
    detray_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-G -src-in-ptx" )
    detray_add_flag( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-lineinfo -src-in-ptx" )

    # Fail on warnings, if asked for that behaviour.
    if(DETRAY_FAIL_ON_WARNINGS)
        if(
            ("${CUDAToolkit_VERSION}" VERSION_GREATER_EQUAL "10.2")
            AND ("${CMAKE_CUDA_COMPILER_ID}" MATCHES "NVIDIA")
        )
            detray_add_flag( CMAKE_CUDA_FLAGS "-Werror all-warnings" )
        elseif("${CMAKE_CUDA_COMPILER_ID}" MATCHES "Clang")
            detray_add_flag( CMAKE_CUDA_FLAGS "-Werror" )
        endif()
    endif()
endif()
