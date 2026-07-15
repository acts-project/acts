# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare tests related to the example code.
traccc_add_test( examples
   "test_options.cpp"
   LINK_LIBRARIES GTest::gtest_main traccc_tests_common traccc::options )
