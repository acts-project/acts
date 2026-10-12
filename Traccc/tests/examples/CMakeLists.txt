# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Declare tests related to the example code.
traccc_add_test( examples
   "test_options.cpp"
   LINK_LIBRARIES GTest::gtest_main traccc_tests_common traccc::options
)
