# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the covfie tests.
detray_add_unit_test( covfie
   "constant_field.cpp"
   LINK_LIBRARIES GTest::gtest_main covfie::core
)
