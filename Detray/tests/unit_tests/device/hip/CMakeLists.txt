# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray HIP unit tests")

cmake_minimum_required(VERSION 3.21) # HIP language support requires minimum 3.21

find_package(HIPToolkit)

# Enable HIP as a language.
enable_language(HIP)

include(detray-compiler-options-hip)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
    detray_add_unit_test(hip_${CMAKE_HIP_PLATFORM}_${algebra}
       "detector_hip_kernel.hpp"
       "detector_hip.cpp"
       "detector_hip_kernel.hip"

       LINK_LIBRARIES GTest::gtest_main HIP::hiprt vecmem::hip detray::core
       detray::algebra_${algebra} detray::test_common detray::test_utils
    )

    # Make hipcc interpret .cpp files as .hip files to make linking work
    # (only needed for the amd backend)
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    )
        set_source_files_properties(detector_hip.cpp PROPERTIES LANGUAGE HIP)
    endif()

    # Add definitions for HIP + algebra
    target_compile_definitions(
        detray_unit_test_hip_${CMAKE_HIP_PLATFORM}_${algebra}
        PRIVATE ${algebra}=${algebra}
    )
endforeach()
