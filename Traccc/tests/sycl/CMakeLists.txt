# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

include(traccc-compiler-options-sycl)

enable_language(SYCL)

traccc_add_test(
    sycl

    test_kalman_fitter_telescope.cpp
    test_ckf_combinatorics_telescope.cpp
    test_ckf_toy_detector.cpp
    test_clusterization.cpp
    test_dpl.sycl
    test_spacepoint_formation.cpp
    test_barrier.sycl
    test_mutex.sycl
    test_unique_lock.sycl
    test_cca.cpp
    test_sanity_contiguous_on.sycl
    test_sanity_ordered_on.sycl
    test_sort.sycl

    LINK_LIBRARIES
    GTest::gtest_main
    vecmem::sycl
    detray::core
    detray::io
    detray::test_common
    traccc::core
    traccc::device_common
    traccc::sycl
    traccc::performance
    traccc::io
    traccc::simulation
    traccc_tests_common
    oneDPL
)
