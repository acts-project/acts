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
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/unmasked.hpp"
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

// Three-dimensional definitions
using test_algebra = test::algebra;
using vector3 = test::vector3;
using point3 = test::point3;
using scalar = test::scalar;
using transform3 = test::transform3;

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

// This defines the local frame test suite
GTEST_TEST(detray_intersection, translated_plane_ray) {
  // Create a shifted plane
  const transform3 shifted(vector3{3.f, 2.f, 10.f});

  // Test ray
  const point3 pos{2.f, 1.f, 0.f};
  const vector3 mom{0.f, 0.f, 1.f};
  const detail::ray<test_algebra> r(pos, 0.f, mom, 0.f);

  // The same test but bound to local frame
  ray_intersector<unmasked<2>, test_algebra, true> pi;
  mask<unmasked<2>, test_algebra> unmasked_bound{};
  const auto hit_bound = pi(r, surface_descriptor<>{}, unmasked_bound, shifted);

  ASSERT_TRUE(hit_bound.is_inside());
  // Global intersection information - unchanged
  const auto global0 =
      unmasked_bound.to_global_frame(shifted, hit_bound.local());
  ASSERT_NEAR(global0[0], 2.f, tol);
  ASSERT_NEAR(global0[1], 1.f, tol);
  ASSERT_NEAR(global0[2], 10.f, tol);
  // Local intersection information
  ASSERT_NEAR(hit_bound.local()[0], -1.f, tol);
  ASSERT_NEAR(hit_bound.local()[1], -1.f, tol);

  // The same test but bound to local frame & masked - inside
  mask<rectangle2D, test_algebra> rect_for_inside{0u, 3.f, 3.f};
  const auto hit_bound_inside =
      pi(r, surface_descriptor<>{}, rect_for_inside, shifted, tol);

  ASSERT_TRUE(hit_bound_inside.is_inside());
  // Global intersection information - unchanged
  const auto global1 =
      rect_for_inside.to_global_frame(shifted, hit_bound_inside.local());
  ASSERT_NEAR(global1[0], 2.f, tol);
  ASSERT_NEAR(global1[1], 1.f, tol);
  ASSERT_NEAR(global1[2], 10.f, tol);
  // Local intersection infoimation - unchanged
  ASSERT_NEAR(hit_bound_inside.local()[0], -1.f, tol);
  ASSERT_NEAR(hit_bound_inside.local()[1], -1.f, tol);

  // The same test but bound to local frame & masked - outside
  mask<rectangle2D, test_algebra> rect_for_outside{0u, 0.5f, 3.5f};
  const auto hit_bound_outside =
      pi(r, surface_descriptor<>{}, rect_for_outside, shifted, tol);
  ASSERT_FALSE(hit_bound_outside.is_inside());
  // Global intersection information - not written out anymore
  const auto global2 =
      rect_for_outside.to_global_frame(shifted, hit_bound_outside.local());
  ASSERT_NEAR(global2[0], 2.f, tol);
  ASSERT_NEAR(global2[1], 1.f, tol);
  ASSERT_NEAR(global2[2], 10.f, tol);
  // Local intersection infoimation - unchanged
  ASSERT_NEAR(hit_bound_outside.local()[0], -1.f, tol);
  ASSERT_NEAR(hit_bound_outside.local()[1], -1.f, tol);
}
