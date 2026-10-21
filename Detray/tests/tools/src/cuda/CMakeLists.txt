# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

message(STATUS "Building detray CUDA command line tools")

# Enable CUDA as a language.
enable_language(CUDA)

include(CMakeFindDependencyMacro)

find_dependency(Boost COMPONENTS program_options REQUIRED)

if(DETRAY_BUILD_TESTING)
    # Build the navigation validation executable.
    detray_add_executable(navigation_validation_cuda
                        "navigation_validation_cuda.cpp"
                        LINK_LIBRARIES GTest::gtest GTest::gtest_main
                        Boost::program_options detray::test_cuda detray::tools
    )

    # Build the material cuda validation executable.
    detray_add_executable(material_validation_cuda
                        "material_validation_cuda.cpp"
                        LINK_LIBRARIES GTest::gtest GTest::gtest_main
                        Boost::program_options detray::test_cuda detray::tools
    )
endif()

if(DETRAY_BUILD_BENCHMARKS)
    # Build benchmarks for multiple algebra plugins
    # Currently vc and smatrix is not supported on device
    set(algebra_plugins "array")
    if(DETRAY_EIGEN_PLUGIN)
        list(APPEND algebra_plugins "eigen")
    endif()

    foreach(algebra ${algebra_plugins})
        detray_add_executable(propagation_benchmark_cuda_${algebra}
        "propagation_benchmark_cuda.cpp"
        LINK_LIBRARIES detray::benchmark_cuda_${algebra} vecmem::cuda detray::tools detray::test_common
        )
    endforeach()
endif()
