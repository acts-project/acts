// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/geometry/shapes/unbounded.hpp"

#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <cassert>
#include <type_traits>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using transform3 = test::transform3;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of an unbounded rectangle shape
GTEST_TEST(detray_masks, unbounded) {
  using shape_t = rectangle2D;
  using unbounded_t = unbounded<shape_t>;

  static_assert(concepts::shape<unbounded_t, test_algebra>);
  static_assert(concepts::planar_shape<unbounded_t, test_algebra>);

  constexpr scalar h{20.f * unit<scalar>::mm};

  mask<unbounded_t, test_algebra> u{0u, h, h};

  // Test local typedefs
  static_assert(std::is_same_v<unbounded_t::shape, shape_t>, "incorrect shape");
  static_assert(std::is_same_v<unbounded_t::boundaries, shape_t::boundaries>,
                "incorrect boundaries");
  static_assert(
      std::is_same_v<unbounded_t::template local_frame_type<test_algebra>,
                     shape_t::template local_frame_type<test_algebra>>,
      "incorrect local frame");

  // Test static members
  EXPECT_TRUE(std::string(unbounded_t::name) ==
              std::string("unbounded rectangle2D"));

  // Test boundary check
  typename mask<unbounded_t, test_algebra>::point3_type p2 = {0.5f, -9.f, 0.f};
  ASSERT_TRUE(u.is_inside(p2, 0.f));

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = u.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(h + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(h + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (h + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (h + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], envelope, tol);

  // Check area
  const scalar a{u.area()};
  EXPECT_NEAR(a, std::numeric_limits<scalar>::max(), tol);
  ASSERT_EQ(a, u.measure());

  const auto centroid = u.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}
