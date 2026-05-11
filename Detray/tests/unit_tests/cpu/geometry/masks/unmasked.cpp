// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/geometry/shapes/unmasked.hpp"

#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of an unmasked plane
GTEST_TEST(detray_masks, unmasked) {
  static_assert(concepts::shape<unmasked<>, test_algebra>);

  point3 p2 = {0.5f, -9.f, 0.f};

  mask<unmasked<>, test_algebra> u{};

  ASSERT_TRUE(u.is_inside(p2, 0.f));

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = u.local_min_bounds(envelope);
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_min_x]));
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_min_y]));
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_min_z]));
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_max_x]));
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_max_y]));
  ASSERT_TRUE(detail::is_invalid_value(loc_bounds[cuboid3D::e_max_z]));

  // Check area
  const scalar a{u.area()};
  EXPECT_NEAR(a, std::numeric_limits<scalar>::max(), tol);
  ASSERT_EQ(a, u.measure());

  const auto centroid = u.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}
