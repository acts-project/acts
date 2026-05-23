// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/utils/find_bound.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// Test upper bound function
GTEST_TEST(detray_utils, upper_bound) {
  std::vector<float> vec = {2.f, 3.f, 5.f, 8.f, 8.f, 8.f, 9.f, 12.f, 12.f};

  auto pos = detray::detail::upper_bound(vec.begin(), vec.end(), 8.f);

  ASSERT_EQ(pos - vec.begin(), 6);
  ASSERT_EQ(*pos, 9.f);

  pos = detray::detail::upper_bound(vec.begin(), vec.end(), 3.f);

  ASSERT_EQ(pos - vec.begin(), 2);
  ASSERT_EQ(*pos, 5.f);

  pos = detray::detail::upper_bound(vec.begin(), vec.end(), 10.f);

  ASSERT_EQ(pos - vec.begin(), 7);
  ASSERT_EQ(*pos, 12.f);
}

// Test lower bound function
GTEST_TEST(detray_utils, lower_bound) {
  std::vector<float> vec = {2.f, 3.f, 5.f, 8.f, 8.f, 8.f, 9.f, 12.f, 12.f};

  auto pos = detray::detail::lower_bound(vec.begin(), vec.end(), 8.f);

  ASSERT_EQ(pos - vec.begin(), 3);
  ASSERT_EQ(*pos, 8.f);

  pos = detray::detail::lower_bound(vec.begin(), vec.end(), 3.f);

  ASSERT_EQ(pos - vec.begin(), 1);
  ASSERT_EQ(*pos, 3.f);

  pos = detray::detail::lower_bound(vec.begin(), vec.end(), 10.f);

  ASSERT_EQ(pos - vec.begin(), 7);
  ASSERT_EQ(*pos, 12.f);
}
