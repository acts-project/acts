// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// detray core
#include "detray/utils/quadratic_equation.hpp"

#include "detray/definitions/algebra.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using scalar = test::scalar;

// This tests the convenience quadratic equation class
GTEST_TEST(detray_utils, quadratic_equation) {
  static constexpr scalar epsilon{1e-5f};

  // No solution
  detail::quadratic_equation<scalar> qe1{1.5f, 0.f, 1.f};

  ASSERT_EQ(qe1.solutions(), 0);

  // One solution
  detail::quadratic_equation<scalar> qe2{1.f, 0.f, 0.f};

  ASSERT_EQ(qe2.solutions(), 1);
  EXPECT_NEAR(qe2.smaller(), 0.f, epsilon);

  detail::quadratic_equation<scalar> qe3{0.f, 1.f, 2.f};

  ASSERT_EQ(qe3.solutions(), 1);
  EXPECT_NEAR(qe3.smaller(), -2.f, epsilon);

  // Two solutions
  detail::quadratic_equation<scalar> qe4{2.f, 5.f, 3.f};

  ASSERT_EQ(qe4.solutions(), 2);
  EXPECT_NEAR(qe4.smaller(), -1.5f, epsilon);
  EXPECT_NEAR(qe4.larger(), -1.f, epsilon);

  detail::quadratic_equation<scalar> qe5{2.f, 5.f, -3.f};

  ASSERT_EQ(qe5.solutions(), 2);
  EXPECT_NEAR(qe5.smaller(), -3.f, epsilon);
  EXPECT_NEAR(qe5.larger(), 0.5f, epsilon);

  detail::quadratic_equation<scalar> qe6{2.f, -5.f, 3.f};

  ASSERT_EQ(qe6.solutions(), 2);
  EXPECT_NEAR(qe6.smaller(), 1.f, epsilon);
  EXPECT_NEAR(qe6.larger(), 1.5f, epsilon);

  detail::quadratic_equation<scalar> qe7{2.f, -5.f, -3.f};

  ASSERT_EQ(qe7.solutions(), 2);
  EXPECT_NEAR(qe7.smaller(), -0.5f, epsilon);
  EXPECT_NEAR(qe7.larger(), 3.f, epsilon);

  detail::quadratic_equation<scalar> qe8{2.f, -5.f, 0.f};

  ASSERT_EQ(qe8.solutions(), 2);
  EXPECT_NEAR(qe8.smaller(), 0.f, epsilon);
  EXPECT_NEAR(qe8.larger(), 5.f / 2.f, epsilon);
}
