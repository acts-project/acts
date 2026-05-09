// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

// Three-dimensional definitions
using test_algebra = test::algebra;
using transform3 = test::transform3;
using vector3 = test::vector3;
using point3 = test::point3;
using point2 = test::point2;
using scalar = test::scalar;

using cartesian = cartesian2D<test_algebra>;
using intersection_t = intersection2D<surface_descriptor<>, test_algebra,
                                      intersection::contains_pos>;
using line_intersector_type =
    ray_intersector<line_circular, test_algebra, true>;

constexpr scalar tol{1e-5f};

// Test simplest case
GTEST_TEST(detray_intersection, line_intersector_case1) {
  // tf3 with Identity rotation and no translation
  const transform3 tf{};

  // Create a track
  std::vector<free_track_parameters<test_algebra>> trks;
  trks.emplace_back(point3{1.f, -1.f, 0.f}, 0.f, vector3{0.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{-1.f, -1.f, 0.f}, 0.f, vector3{0.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{1.f, -1.f, 2.f}, 0.f, vector3{0.f, 1.f, -1.f}, -1.f);

  // Infinite wire with 10 mm radial cell size
  const mask<line_circular, test_algebra> ln{
      0u, 10.f, std::numeric_limits<scalar>::infinity()};

  // Test intersect
  std::vector<intersection_t> is(3u);
  is[0] = line_intersector_type()(detail::ray(trks[0]), surface_descriptor<>{},
                                  ln, tf, tol);
  is[1] = line_intersector_type()(detail::ray(trks[1]), surface_descriptor<>{},
                                  ln, tf, tol);
  is[2] = line_intersector_type()(detail::ray(trks[2]), surface_descriptor<>{},
                                  ln, tf, tol);

  EXPECT_TRUE(is[0].is_inside());
  EXPECT_EQ(is[0].path(), 1.f);

  const point3 global0 = ln.to_global_frame(tf, is[0].local());
  point3 x{1.f, 0.f, 0.f};
  EXPECT_NEAR(global0[0], x[0], tol);
  EXPECT_NEAR(global0[1], x[1], tol);
  EXPECT_NEAR(global0[2], x[2], tol);
  EXPECT_EQ(is[0].local()[0], -1.f);  // right
  EXPECT_EQ(is[0].local()[1], 0.f);

  EXPECT_TRUE(is[1].is_inside());
  EXPECT_EQ(is[1].path(), 1.f);
  const auto global1 = ln.to_global_frame(tf, is[1].local());
  EXPECT_NEAR(global1[0], -1.f, tol);
  EXPECT_NEAR(global1[1], 0.f, tol);
  EXPECT_NEAR(global1[2], 0.f, tol);
  EXPECT_EQ(is[1].local()[0], 1.f);  // left
  EXPECT_EQ(is[1].local()[1], 0.f);

  EXPECT_TRUE(is[2].is_inside());
  EXPECT_NEAR(is[2].path(), constant<scalar>::sqrt2, tol);
  const auto global2 = ln.to_global_frame(tf, is[2].local());
  EXPECT_NEAR(global2[0], 1.f, tol);
  EXPECT_NEAR(global2[1], 0.f, tol);
  EXPECT_NEAR(global2[2], 1.f, tol);
  EXPECT_NEAR(is[2].local()[0], -1.f, tol);  // right
  EXPECT_NEAR(is[2].local()[1], 1.f, tol);
}

// Test inclined wire
GTEST_TEST(detray_intersection, line_intersector_case2) {
  // tf3 with skewed axis
  const vector3 x{1.f, 0.f, -1.f};
  const vector3 z{1.f, 0.f, 1.f};
  const vector3 t{1.f, 1.f, 1.f};
  const transform3 tf{t, vector::normalize(z), vector::normalize(x)};

  // Create a track
  const point3 pos{1.f, -1.f, 0.f};
  const vector3 dir{0.f, 1.f, 0.f};
  const free_track_parameters<test_algebra> trk(pos, 0.f, dir, -1.f);

  // Infinite wire with 10 mm
  // radial cell size
  const mask<line_circular, test_algebra> ln{
      0u, 10.f, std::numeric_limits<scalar>::infinity()};

  // Test intersect
  const intersection_t is = line_intersector_type()(
      detail::ray<test_algebra>(trk), surface_descriptor<>{}, ln, tf, tol);

  EXPECT_TRUE(is.is_inside());
  EXPECT_NEAR(is.path(), 2.f, tol);
  const auto global = ln.to_global_frame(tf, is.local());
  EXPECT_NEAR(global[0], 1.f, tol);
  EXPECT_NEAR(global[1], 1.f, tol);
  EXPECT_NEAR(global[2], 0.f, tol);
  EXPECT_NEAR(is.local()[0], -constant<scalar>::inv_sqrt2, tol);  // right
  EXPECT_NEAR(is.local()[1], -constant<scalar>::inv_sqrt2, tol);
}

GTEST_TEST(detray_intersection, line_intersector_square_scope) {
  // tf3 with Identity rotation and no translation
  const transform3 tf{};

  /// Create a track
  std::vector<free_track_parameters<test_algebra>> trks;
  trks.emplace_back(point3{2.f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{1.9f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{2.1f, 0.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f}, -1.f);

  trks.emplace_back(point3{-2.f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{-1.9f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{-2.1f, 0.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);

  trks.emplace_back(point3{0.f, -2.f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{0.f, -1.9f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                    -1.f);
  trks.emplace_back(point3{0.f, -2.1f, 0.f}, 0.f, vector3{-1.f, 1.f, 0.f},
                    -1.f);

  trks.emplace_back(point3{0.f, -2.f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{0.f, -1.9f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);
  trks.emplace_back(point3{0.f, -2.1f, 0.f}, 0.f, vector3{1.f, 1.f, 0.f}, -1.f);

  // Infinite wire with 1 mm square cell size
  mask<line_square, test_algebra, std::uint_least16_t> ln{
      0u, 1.f, std::numeric_limits<scalar>::infinity()};

  // Test intersect
  std::vector<intersection_t> is;
  for (const auto& trk : trks) {
    is.push_back(line_intersector_type()(detail::ray<test_algebra>(trk),
                                         surface_descriptor<>{}, ln, tf, tol));
  }

  ASSERT_TRUE(is.size() >= 12);

  EXPECT_TRUE(is[0].is_inside());
  EXPECT_NEAR(is[0].path(), constant<scalar>::sqrt2, tol);
  const auto local0 = ln.to_local_frame3D(
      tf,
      detail::ray(trks[0]).pos() + is[0].path() * detail::ray(trks[0]).dir(),
      detail::ray(trks[0]).dir());
  const auto global0 = ln.to_global_frame(tf, local0);
  EXPECT_NEAR(global0[0], 1.f, tol);
  EXPECT_NEAR(global0[1], 1.f, tol);
  EXPECT_NEAR(global0[2], 0.f, tol);
  EXPECT_NEAR(is[0].local()[0], -constant<scalar>::sqrt2, tol);
  EXPECT_NEAR(is[0].local()[1], 0.f, tol);

  EXPECT_TRUE(is[1].is_inside());
  EXPECT_TRUE(std::signbit(is[1].local()[0]));
  EXPECT_FALSE(is[2].is_inside());
  EXPECT_TRUE(std::signbit(is[2].local()[0]));

  EXPECT_TRUE(is[3].is_inside());
  EXPECT_FALSE(std::signbit(is[3].local()[0]));
  EXPECT_TRUE(is[4].is_inside());
  EXPECT_FALSE(std::signbit(is[4].local()[0]));
  EXPECT_FALSE(is[5].is_inside());
  EXPECT_FALSE(std::signbit(is[5].local()[0]));

  EXPECT_TRUE(is[6].is_inside());
  EXPECT_FALSE(std::signbit(is[6].local()[0]));
  EXPECT_TRUE(is[7].is_inside());
  EXPECT_FALSE(std::signbit(is[7].local()[0]));
  EXPECT_FALSE(is[8].is_inside());
  EXPECT_FALSE(std::signbit(is[8].local()[0]));

  EXPECT_TRUE(is[9].is_inside());
  EXPECT_TRUE(std::signbit(is[9].local()[0]));
  EXPECT_TRUE(is[10].is_inside());
  EXPECT_TRUE(std::signbit(is[10].local()[0]));
  EXPECT_FALSE(is[11].is_inside());
  EXPECT_TRUE(std::signbit(is[11].local()[0]));
}
