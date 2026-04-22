// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/geometry/shapes/trapezoid2D.hpp"

#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/ratio_test.hpp"

// GTest include
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of a trapezoid
GTEST_TEST(detray_masks, trapezoid2D) {
  static_assert(concepts::shape<trapezoid2D, test_algebra>);
  static_assert(concepts::planar_shape<trapezoid2D, test_algebra>);

  point3 p2_in = {1.f, -0.5f, 0.f};
  point3 p2_edge = {2.5f, 1.f, 0.f};
  point3 p2_out = {3.f, 1.5f, 0.f};

  constexpr scalar hx_miny{1.f * unit<scalar>::mm};
  constexpr scalar hx_maxy{3.f * unit<scalar>::mm};
  constexpr scalar hy{2.f * unit<scalar>::mm};
  constexpr scalar divisor{1.f / (2.f * hy)};

  mask<trapezoid2D, test_algebra> t2{0u, hx_miny, hx_maxy, hy, divisor};

  ASSERT_NEAR(t2[trapezoid2D::e_half_length_0], hx_miny, tol);
  ASSERT_NEAR(t2[trapezoid2D::e_half_length_1], hx_maxy, tol);
  ASSERT_NEAR(t2[trapezoid2D::e_half_length_2], hy, tol);
  ASSERT_NEAR(t2[trapezoid2D::e_divisor], divisor, tol);

  ASSERT_TRUE(t2.is_inside(p2_in));
  ASSERT_TRUE(t2.is_inside(p2_edge));
  ASSERT_FALSE(t2.is_inside(p2_out));
  // Move outside point inside using a tolerance

  // Check area
  const scalar a{t2.area()};
  EXPECT_NEAR(a, 16.f * unit<scalar>::mm2, tol);
  ASSERT_EQ(a, t2.measure());

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = t2.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(hx_maxy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (hx_maxy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], envelope, tol);

  const auto centroid = t2.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 1.f / 3.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, trapezoid2D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<trapezoid2D, test_algebra> &tp,
                    const test::transform3 &trf, const scalar t) {
      return tp.is_inside(trf, p, t);
    }
  };

  constexpr mask<trapezoid2D, test_algebra> tp{0u, 2.f, 3.f, 4.f,
                                               1.f / (2.f * 4.f)};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  // x- and y-coordinates yield a valid local position on the underlying plane
  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, tp, trf, t);

  const scalar area{tp.measure()};
  const scalar world{size * size};

  ASSERT_NEAR(ratio, area / world, 0.004f);
}
