# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

enable_language(CUDA)

include(traccc-compiler-options-cuda)

find_package(CUDAToolkit REQUIRED)

add_library(
    traccc_tests_cuda_main
    STATIC
    cuda_main.cpp
)

target_link_libraries(
    traccc_tests_cuda_main
    PRIVATE
    CUDA::cudart
    GTest::gtest
)

# Set up Thrust specifically for the CUDA unit tests.
thrust_create_target( traccc_test_cuda_thrust
  HOST CPP
  DEVICE CUDA )

traccc_add_test(
    cuda

    # Define the sources for the test.
    test_ambiguity_resolution.cpp
    test_basic.cu
    test_barrier.cu
    test_cca.cpp
    test_ckf_combinatorics_telescope.cpp
    test_ckf_toy_detector.cpp
    test_copy.cu
    test_kalman_fitter_telescope.cpp
    test_clusterization.cpp
    test_copy.cu
    test_grid2.cpp
    test_grid2.cu
    test_spacepoint_formation.cpp
    test_thrust.cu
    test_sync.cu
    test_array_wrapper.cu
    test_mutex.cu
    test_unique_lock.cu
    test_sanity_contiguous_on.cu
    test_sanity_ordered_on.cu
    test_sort.cu

    LINK_LIBRARIES
    CUDA::cudart
    traccc_test_cuda_thrust
    GTest::gtest
    vecmem::cuda
    detray::core
    detray::io
    detray::test_common
    traccc::core
    traccc::device_common
    traccc::cuda_utils
    traccc::cuda
    traccc::performance
    traccc::io
    traccc::simulation
    traccc_tests_cuda_main
    traccc_tests_common
)

# For CUDA 11 turn on separable compilation. This is necessary for using
# Thrust 2.1.0.
if( ( "${CMAKE_CUDA_COMPILER_ID}" STREQUAL "NVIDIA" ) AND
    ( "${CMAKE_CUDA_COMPILER_VERSION}" VERSION_LESS "12.0" ) )
    set_target_properties( traccc_test_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON )
endif()
