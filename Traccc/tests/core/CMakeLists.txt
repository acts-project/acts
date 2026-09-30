# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Declare the core library test(s).
traccc_add_test(core
   "test_algorithm.cpp"
   "test_axis.cpp"
   "test_copy.cpp"
   "test_grid2.cpp"
   "test_populator.cpp"
   "test_pvalue.cpp"
   "test_ranges.cpp"
   "test_serializer.cpp"
   "test_simulation.cpp"
   "test_subspace.cpp"
   "particle.cpp"
   LINK_LIBRARIES GTest::gtest_main vecmem::core traccc_tests_common
   traccc::core traccc::io traccc::performance traccc::simulation
)

if(TRACCC_USE_ROOT)
    find_package(ROOT COMPONENTS Core RIO Hist REQUIRED)
    target_link_libraries(traccc_test_core PRIVATE ROOT::Core)
    target_compile_definitions(traccc_test_core PRIVATE TRACCC_HAVE_ROOT)
endif()
