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
#include "detray/geometry/shapes/rectangle2D.hpp"

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

constexpr scalar hx{1.f * unit<scalar>::mm};
constexpr scalar hy{9.3f * unit<scalar>::mm};
constexpr scalar hz{0.5f * unit<scalar>::mm};

/// This tests the basic functionality of a rectangle
GTEST_TEST(detray_masks, rectangle2D) {
  static_assert(concepts::shape<rectangle2D, test_algebra>);
  static_assert(concepts::rectilinear_shape<rectangle2D, test_algebra>);
  static_assert(concepts::planar_shape<rectangle2D, test_algebra>);

  point3 p2_in = {0.5f, -9.f, 0.f};
  point3 p2_edge = {1.f, 9.3f, 0.f};
  point3 p2_out = {1.5f, -9.f, 0.f};

  mask<rectangle2D, test_algebra> r2{0u, hx, hy};

  ASSERT_NEAR(r2[rectangle2D::e_half_x], hx, tol);
  ASSERT_NEAR(r2[rectangle2D::e_half_y], hy, tol);

  ASSERT_TRUE(r2.is_inside(p2_in));
  ASSERT_TRUE(r2.is_inside(p2_edge));
  ASSERT_FALSE(r2.is_inside(p2_out));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(r2.is_inside(p2_out, 1.f));

  // Check area
  const scalar a{r2.area()};
  EXPECT_NEAR(a, 37.2f * unit<scalar>::mm2, tol);
  ASSERT_EQ(a, r2.measure());

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = r2.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(hx + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (hx + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], envelope, tol);

  const auto centroid = r2.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, rectangle2D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<rectangle2D, test_algebra> &r,
                    const test::transform3 &trf, const scalar t) {
      return r.is_inside(trf, p, t);
    }
  };

  constexpr mask<rectangle2D, test_algebra> r{0u, 3.f, 4.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  // x- and y-coordinates yield a valid local position on the underlying plane
  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, r, trf, t);

  const scalar area{r.measure()};
  const scalar world{size * size};

  ASSERT_NEAR(ratio, area / world, 0.02f);
}

/// This tests the basic functionality of a cuboid3D
GTEST_TEST(detray_masks, cuboid3D) {
  static_assert(concepts::shape<cuboid3D, test_algebra>);
  static_assert(concepts::rectilinear_shape<cuboid3D, test_algebra>);

  point3 p2_in = {0.5f, 8.0f, -0.4f};
  point3 p2_edge = {1.f, 9.3f, 0.5f};
  point3 p2_out = {1.5f, -9.f, 0.55f};

  mask<cuboid3D, test_algebra> c3{0u, -hx, -hy, -hz, hx, hy, hz};

  ASSERT_NEAR(c3[cuboid3D::e_min_x], -hx, tol);
  ASSERT_NEAR(c3[cuboid3D::e_min_y], -hy, tol);
  ASSERT_NEAR(c3[cuboid3D::e_min_z], -hz, tol);
  ASSERT_NEAR(c3[cuboid3D::e_max_x], hx, tol);
  ASSERT_NEAR(c3[cuboid3D::e_max_y], hy, tol);
  ASSERT_NEAR(c3[cuboid3D::e_max_z], hz, tol);

  ASSERT_TRUE(c3.is_inside(p2_in));
  ASSERT_TRUE(c3.is_inside(p2_edge));
  ASSERT_FALSE(c3.is_inside(p2_out));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(c3.is_inside(p2_out, 1.f));

  // Check volume
  const scalar v{c3.volume()};
  EXPECT_NEAR(v, 37.2f * unit<scalar>::mm3, tol);
  ASSERT_EQ(v, c3.measure());

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = c3.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(hx + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -(hz + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (hx + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (hy + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], (hz + envelope), tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, cuboid3D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<cuboid3D, test_algebra> &cb,
                    const test::transform3 &trf, const scalar t) {
      const point3 loc_p{cb.to_local_frame3D(trf, p)};
      return cb.is_inside(loc_p, t);
    }
  };

  constexpr mask<cuboid3D, test_algebra> cb{0u, 0.f, 0.f, 0.f, 3.f, 4.f, 1.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, cb, trf, t);

  const scalar volume{cb.measure()};
  const scalar world{size * size * size};

  ASSERT_NEAR(ratio, volume / world, 0.002f);
}
