# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include(traccc-compiler-options-sycl)

include(vecmem-check-language)
vecmem_check_language(SYCL)
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
