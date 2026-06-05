# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

enable_language(HIP)
traccc_add_test(
    hip
    # Define the sources for the test.
    test_thrust.hip
    LINK_LIBRARIES
    roc::rocthrust
    GTest::gtest_main
    vecmem::core
    vecmem::hip
)

set_target_properties( traccc_test_hip PROPERTIES
   POSITION_INDEPENDENT_CODE TRUE )
