# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Helper function for running a test in a specific directory.
function(_run_test_in_dir test dir)
    # Make the directory.
    file(MAKE_DIRECTORY "${dir}")
    # Run the test in that directory.
    set_tests_properties(
        "detray_integration_test_${test}"
        PROPERTIES WORKING_DIRECTORY "${dir}"
    )
endfunction()

detray_add_integration_test( io_roundtrip
    "io_json_detector_roundtrip.cpp"
    LINK_LIBRARIES GTest::gtest_main vecmem::core detray::core_array
    detray::io_array detray::test_common detray::test_cpu
)

# Run integration tests only after unit tests passed
set_tests_properties(
    detray_integration_test_io_roundtrip
    PROPERTIES DEPENDS "detray_unit_test_io_payloads;detray_unit_test_io_writer"
)

_run_test_in_dir( io_roundtrip
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/io_roundtrip_test_rundir"
)
