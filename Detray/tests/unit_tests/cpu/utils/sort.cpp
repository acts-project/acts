// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/sort.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// Test sort functions
GTEST_TEST(detray_utils, insertion_sort) {
  std::vector<double> vec = {4.1, 5., 1.2, 1.4, 9.};
  std::vector<double> vec_sorted = {1.2, 1.4, 4.1, 5., 9.};

  detray::detail::insertion_sort(vec.begin(), vec.end());

  ASSERT_EQ(vec, vec_sorted);
}

GTEST_TEST(detray_utils, selection_sort) {
  std::vector<double> vec = {4.1, 5., 1.2, 1.4, 9.};
  std::vector<double> vec_sorted = {1.2, 1.4, 4.1, 5., 9.};

  detray::detail::selection_sort(vec.begin(), vec.end());

  ASSERT_EQ(vec, vec_sorted);
}
