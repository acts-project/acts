# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Declare the io library test(s).
traccc_add_test( io
   "test_bfield.cpp"
   "test_csv.cpp"
   "test_event_data.cpp"
   "test_json.cpp"
   LINK_LIBRARIES GTest::gtest_main traccc_tests_common
                  traccc::core traccc::io traccc::performance
)

target_compile_definitions(
    traccc_test_io
    PRIVATE TRACCC_TEST_IO_MOCK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/mock_data"
)
