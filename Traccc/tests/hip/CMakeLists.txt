# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2024-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable the HIP language.
enable_language(HIP)

add_library(traccc_test_hip_device STATIC test_thrust.hip)
target_link_libraries(
    traccc_test_hip_device
    PRIVATE roc::rocthrust GTest::gtest vecmem::core vecmem::hip
)
set_target_properties(
    traccc_test_hip_device
    PROPERTIES POSITION_INDEPENDENT_CODE TRUE
)

traccc_add_test(
    hip
    # Define the sources for the test.
    test_cca.cpp
    LINK_LIBRARIES
    GTest::gtest_main
    vecmem::core
    vecmem::hip
    traccc::core
    traccc::device_common
    traccc::hip_utils
    traccc::hip
    traccc_tests_common
    traccc_test_hip_device
)

set_target_properties(traccc_test_hip PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
