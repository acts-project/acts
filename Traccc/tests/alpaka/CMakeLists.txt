# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include(traccc-alpaka-functions)
traccc_enable_language_alpaka()

traccc_add_alpaka_test( alpaka
   alpaka_basic.cpp
   LINK_LIBRARIES
   GTest::gtest_main
   traccc_tests_common
   alpaka::alpaka
   vecmem::core
   traccc::alpaka
)

traccc_add_test(alpaka_cxx
   test_cca.cpp
   test_ckf_combinatorics_telescope.cpp
   test_ckf_toy_detector.cpp
   LINK_LIBRARIES
   GTest::gtest_main
   traccc_tests_common
   vecmem::core
   traccc::alpaka
   traccc::io
   traccc::performance
   traccc::simulation
)
