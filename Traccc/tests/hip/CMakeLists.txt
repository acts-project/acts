# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable the HIP language.
enable_language(HIP)

find_package(HIPToolkit REQUIRED)

traccc_add_test(
    hip
    # Define the sources for the test.
    test_thrust.hip
    test_cca.cpp
    LINK_LIBRARIES
    HIP::hiprt
    roc::rocthrust
    GTest::gtest_main
    vecmem::core
    vecmem::hip
    traccc::core
    traccc::device_common
    traccc::hip
    traccc_tests_common
)

set_target_properties( traccc_test_hip PROPERTIES
   POSITION_INDEPENDENT_CODE TRUE )
