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
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

namespace {

// Three-dimensional definitions
using test_algebra = test::algebra;
using transform3_t = test::transform3;
using vector3 = test::vector3;
using point3 = test::point3;
using scalar = test::scalar;
using ray_t = detray::detail::ray<test_algebra>;

constexpr scalar not_defined = std::numeric_limits<scalar>::max();
constexpr scalar tol{1e-5f};

const scalar r{4.f};
const scalar hz{10.f};

}  // anonymous namespace

// This checks both solutions of a ray-cylinder intersection
GTEST_TEST(detray_intersection, translated_cylinder) {
  // Create a translated cylinder and test untersection
  const transform3_t shifted(vector3{3.f, 2.f, 10.f});
  ray_intersector<cylinder2D, test_algebra, true> ci;

  // Test ray
  const point3 ori = {3.f, 2.f, 5.f};
  const point3 dir = {1.f, 0.f, 0.f};
  ray_t ray(ori, 0.f, dir, 0.f);

  // Intersect:
  mask<cylinder2D, test_algebra, std::uint_least16_t> cylinder{0u, r, -hz, hz};
  const auto hits_bound =
      ci(ray, surface_descriptor<>{}, cylinder, shifted, tol, -not_defined);

  // first intersection lies behind the track
  EXPECT_TRUE(hits_bound[0].is_inside());
  ASSERT_FALSE(hits_bound[0].is_along());

  const auto global0 = cylinder.to_global_frame(shifted, hits_bound[0].local());
  EXPECT_NEAR(global0[0], -1.f, tol);
  EXPECT_NEAR(global0[1], 2.f, tol);
  EXPECT_NEAR(global0[2], 5.f, tol);

  ASSERT_TRUE(hits_bound[0].local()[0] != not_defined &&
              hits_bound[0].local()[1] != not_defined);
  // p2[0] = r * phi : 180deg in the opposite direction with r = 4
  EXPECT_NEAR(hits_bound[0].local()[0], 4.f * constant<scalar>::pi, tol);
  EXPECT_NEAR(hits_bound[0].local()[1], -5., tol);

  // second intersection lies in front of the track
  EXPECT_TRUE(hits_bound[1].is_inside());
  EXPECT_TRUE(hits_bound[1].is_along());

  const auto global1 = cylinder.to_global_frame(shifted, hits_bound[1].local());
  EXPECT_NEAR(global1[0], 7.f, tol);
  EXPECT_NEAR(global1[1], 2.f, tol);
  EXPECT_NEAR(global1[2], 5.f, tol);

  ASSERT_TRUE(hits_bound[1].local()[0] != not_defined &&
              hits_bound[1].local()[1] != not_defined);
  EXPECT_NEAR(hits_bound[1].local()[0], 0.f, tol);
  EXPECT_NEAR(hits_bound[1].local()[1], -5.f, tol);
}

// This checks the solution of a ray-cylinder portal intersection against
// those obtained from the general cylinder intersector.
GTEST_TEST(detray_intersection, cylinder_portal) {
  // Test ray
  const point3 ori = {1.f, 0.5f, 1.f};
  const point3 dir = vector::normalize(vector3{1.f, 1.f, 1.f});
  const ray_t ray(ori, 0.f, dir, 0.f);

  // Create a concentric cylinder and test intersection
  const transform3_t identity{};
  mask<concentric_cylinder2D, test_algebra, std::uint_least16_t> cylinder{
      0u, r, -hz, hz};

  ray_intersector<cylinder2D, test_algebra, true> ci;
  ray_intersector<concentric_cylinder2D, test_algebra, true> cpi;

  // Intersect
  const auto hits_cylindrical =
      ci(ray, surface_descriptor<>{}, cylinder, identity, tol);
  const auto hit_cocylindrical =
      cpi(ray, surface_descriptor<>{}, cylinder, identity, tol);

  ASSERT_TRUE(hits_cylindrical[1].is_inside());
  ASSERT_TRUE(hit_cocylindrical.is_inside());
  ASSERT_TRUE(hits_cylindrical[1].is_along());
  ASSERT_TRUE(hit_cocylindrical.is_along());

  const auto global0 =
      cylinder.to_global_frame(identity, hits_cylindrical[1].local());
  const auto global1 =
      cylinder.to_global_frame(identity, hit_cocylindrical.local());
  EXPECT_NEAR(vector::perp(global0), r, tol);
  EXPECT_NEAR(vector::perp(global1), r, tol);

  EXPECT_NEAR(global0[0], global1[0], tol);
  EXPECT_NEAR(global0[1], global1[1], tol);
  EXPECT_NEAR(global0[2], global1[2], tol);
  ASSERT_TRUE(hits_cylindrical[1].local()[0] != not_defined &&
              hits_cylindrical[1].local()[1] != not_defined);
  ASSERT_TRUE(hit_cocylindrical.local()[0] != not_defined &&
              hit_cocylindrical.local()[1] != not_defined);
  EXPECT_NEAR(hits_cylindrical[1].local()[0], hit_cocylindrical.local()[0],
              tol);
  EXPECT_NEAR(hits_cylindrical[1].local()[1], hit_cocylindrical.local()[1],
              tol);
}
