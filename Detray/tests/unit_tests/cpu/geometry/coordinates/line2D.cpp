// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/geometry/coordinates/line2D.hpp"

#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point2 = test::point2;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;

constexpr scalar isclose{1e-5f};

GTEST_TEST(detray_coordinates, line2D_case1) {
  // Preparation work
  vector3 z = {1.f, 1.f, 1.f};
  z = vector::normalize(z);
  vector3 x = {1.f, 0.f, -1.f};
  x = vector::normalize(x);
  const point3 t = {0.f, 0.f, 0.f};
  const transform3 trf(t, z, x);
  const point3 global1 = {1.f, 1.5f, 0.5f};
  const vector3 mom = {0.f, 1.f, 1.f};
  const vector3 d = vector::normalize(mom);

  const line2D<test_algebra> l2;

  static_assert(concepts::coordinate_frame<line2D<test_algebra>>);
  static_assert(concepts::line_frame<line2D<test_algebra>>);

  // Global to local transformation
  const point3 local = l2.global_to_local_3D(trf, global1, d);

  // Check if the local position is correct
  ASSERT_NEAR(local[0], -constant<scalar>::inv_sqrt2, isclose);
  ASSERT_NEAR(local[1], std::sqrt(3.f), isclose);

  // Local to global transformation
  const point3 global2 = l2.local_to_global(trf, local);

  // Check if the same global position is obtained
  ASSERT_NEAR(global1[0], global2[0], isclose);
  ASSERT_NEAR(global1[1], global2[1], isclose);
  ASSERT_NEAR(global1[2], global2[2], isclose);

  // Normal vector
  const vector3 n = l2.normal(trf);
  ASSERT_NEAR(n[0], z[0], isclose);
  ASSERT_NEAR(n[1], z[1], isclose);
  ASSERT_NEAR(n[2], z[2], isclose);
}

GTEST_TEST(detray_coordinates, line2D_case2) {
  // Preparation work
  vector3 z = {1.f, 2.f, 3.f};
  z = vector::normalize(z);
  vector3 x = {2.f, -4.f, 2.f};
  x = vector::normalize(x);
  const point3 t = {0.f, 0.f, 0.f};
  const transform3 trf(t, z, x);
  const point2 local1 = {1.f, 2.f};
  const vector3 mom = {1.f, 6.f, -2.f};
  const vector3 d = vector::normalize(mom);
  struct dummy_mask {
  } mask;

  const line2D<test_algebra> l2;

  static_assert(concepts::coordinate_frame<line2D<test_algebra>>);
  static_assert(concepts::line_frame<line2D<test_algebra>>);

  // local to global transformation
  const point3 global = l2.local_to_global(trf, mask, local1, d);

  // global to local transform
  const point3 local2 = l2.global_to_local_3D(trf, global, d);

  // Check if the same local position is obtained
  ASSERT_NEAR(local1[0], local2[0], isclose);
  ASSERT_NEAR(local1[1], local2[1], isclose);
}
