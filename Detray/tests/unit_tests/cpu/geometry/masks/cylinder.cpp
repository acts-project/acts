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
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder3D.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/ratio_test.hpp"

// GTest include
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

constexpr scalar tol{1e-4f};

constexpr scalar r{3.f * unit<scalar>::mm};
constexpr scalar hz{4.f * unit<scalar>::mm};

/// This tests the basic functionality of a 2D cylinder
GTEST_TEST(detray_masks, cylinder2D) {
  static_assert(concepts::shape<cylinder2D, test_algebra>);
  static_assert(concepts::cylindrical_shape<cylinder2D, test_algebra>);

  point3 p2_in = {r, static_cast<scalar>(-1.f), r};
  point3 p2_edge = {r, hz, r};
  point3 p2_out = {3.5f, 4.5f, 3.5f};

  mask<cylinder2D, test_algebra> c{0u, r, -hz, hz};

  ASSERT_NEAR(c[cylinder2D::e_r], r, tol);
  ASSERT_NEAR(c[cylinder2D::e_lower_z], -hz, tol);
  ASSERT_NEAR(c[cylinder2D::e_upper_z], hz, tol);

  ASSERT_TRUE(c.is_inside(p2_in));
  ASSERT_TRUE(c.is_inside(p2_edge));
  ASSERT_FALSE(c.is_inside(p2_out));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(c.is_inside(p2_out, 0.6f));

  // Check area
  const scalar a{c.area()};
  EXPECT_NEAR(a, 150.796447f * unit<scalar>::mm2, tol);
  ASSERT_EQ(a, c.measure());

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = c.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -(hz + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], (hz + envelope), tol);

  const auto centroid = c.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, cylinder2D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<cylinder2D, test_algebra> &cyl,
                    const test::transform3 &trf, const scalar t) {
      return cyl.is_inside(trf, p, t);
    }
  };

  constexpr mask<cylinder2D, test_algebra> cyl{0u, r, 0.f, 5.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(10000, 2))};

  // Make sure the test point is on the cylinder before calling the
  // mask(this is normally ensured by the intersector)
  std::vector<point3> points = test::generate_regular_points<cylinder2D>(
      n_points, {cyl[cylinder2D::e_r], size});

  scalar ratio = test::ratio_test<mask_check>(points, cyl, trf, t);

  const scalar area{cyl.measure()};
  const scalar world{2.f * constant<scalar>::pi * r * size};

  ASSERT_NEAR(ratio, area / world, 0.02f);
}

/// This tests the basic functionality of a 3D cylinder
GTEST_TEST(detray_masks, cylinder3D) {
  static_assert(concepts::shape<cylinder3D, test_algebra>);
  static_assert(concepts::cylindrical_shape<cylinder3D, test_algebra>);

  point3 p3_in = {r, static_cast<scalar>(0.f), static_cast<scalar>(-1.f)};
  point3 p3_edge = {static_cast<scalar>(0.f), r, hz};
  point3 p3_out = {r * constant<scalar>::inv_sqrt2,
                   r * constant<scalar>::inv_sqrt2, static_cast<scalar>(4.5f)};
  point3 p3_off = {1.f, 1.f, -9.f};

  mask<cylinder3D, test_algebra> c{
      0u, 0.f, -constant<scalar>::pi, -hz, r, constant<scalar>::pi, hz};

  ASSERT_NEAR(c[cylinder3D::e_min_r], 0.f, tol);
  ASSERT_NEAR(c[cylinder3D::e_max_r], r, tol);
  ASSERT_NEAR(c[cylinder3D::e_min_phi], -constant<scalar>::pi, tol);
  ASSERT_NEAR(c[cylinder3D::e_max_phi], constant<scalar>::pi, tol);
  ASSERT_NEAR(c[cylinder3D::e_min_z], -hz, tol);
  ASSERT_NEAR(c[cylinder3D::e_max_z], hz, tol);

  ASSERT_TRUE(c.is_inside(p3_in));
  ASSERT_TRUE(c.is_inside(p3_edge));
  ASSERT_FALSE(c.is_inside(p3_out));
  ASSERT_FALSE(c.is_inside(p3_off));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(c.is_inside(p3_out, 0.6f));

  // Check volume
  const scalar v{c.volume()};
  EXPECT_NEAR(v, 226.194671058f * unit<scalar>::mm3, tol);
  ASSERT_EQ(v, c.measure());

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = c.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -(hz + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (r + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], (hz + envelope), tol);

  const auto centroid = c.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, cylinder3D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<cylinder3D, test_algebra> &cyl,
                    const test::transform3 &trf, const scalar t) {
      return cyl.is_inside(trf, p, t);
    }
  };

  constexpr mask<cylinder3D, test_algebra> cyl{
      0u, 1.f, -constant<scalar>::pi, 0.f, 3.f, constant<scalar>::pi, 5.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, cyl, trf, t);

  const scalar volume{cyl.measure()};
  const scalar world{size * size * size};

  ASSERT_NEAR(ratio, volume / world, 0.001f);
}
