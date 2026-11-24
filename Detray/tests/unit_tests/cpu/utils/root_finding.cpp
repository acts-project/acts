// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/utils/root_finding.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <utility>

using namespace detray;

using scalar = test::scalar;

namespace {

// f(s) = (s + 2) * ((s - 11)^2 + 1): single root at s = -2, f(0) = 244
auto evaluate_func = [](const scalar s) {
  const scalar q{(s - 11.f) * (s - 11.f) + 1.f};
  return std::make_pair((s + 2.f) * q, q + 2.f * (s + 2.f) * (s - 11.f));
};

constexpr scalar tol{1e-3f};

}  // namespace

// A start value of exactly 0 must not be returned as root unchecked
GTEST_TEST(detray_utils, newton_raphson_start_at_zero) {
  auto func = evaluate_func;
  const auto [s, err] = newton_raphson(func, scalar{0.f});

  EXPECT_NEAR(s, -2.f, tol);
}

// Starting at s = 10, the bracket expands from [8, 12] to [-16, 16], so the
// iteration starts at its midpoint s = 0
GTEST_TEST(detray_utils, newton_raphson_safe_bracket_midpoint_at_zero) {
  auto func = evaluate_func;

  darray<scalar, 2> br{};
  auto f = [&func](const scalar x) { return func(x).first; };
  ASSERT_TRUE(expand_bracket(scalar{8.f}, scalar{12.f}, f, br));
  ASSERT_EQ(br[0], -16.f);
  ASSERT_EQ(br[1], 16.f);

  const auto [s, err] = newton_raphson_safe(func, scalar{10.f});

  EXPECT_NEAR(s, -2.f, tol);
}
