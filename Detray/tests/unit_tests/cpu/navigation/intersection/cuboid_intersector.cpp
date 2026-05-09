// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/navigation/intersection/bounding_box/cuboid_intersector.hpp"

#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/bounding_volume.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include
#include <gtest/gtest.h>

using namespace detray;

namespace {

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;
using point3 = test::point3;

// cuboid
constexpr scalar x_min{1.f * unit<scalar>::mm};
constexpr scalar x_max{3.f * unit<scalar>::mm};
constexpr scalar y_min{0.f * unit<scalar>::mm};
constexpr scalar y_max{2.f * unit<scalar>::mm};
constexpr scalar z_min{2.f * unit<scalar>::mm};
constexpr scalar z_max{3.f * unit<scalar>::mm};

// envelope around wrapped object (scalor in percent)
constexpr scalar envelope{1.01f};

}  // anonymous namespace

// This test the intersection between ray and cuboid aabb
GTEST_TEST(detray_intersection, cuboid_aabb_intersector) {
  // Test ray
  const point3 pos{2.f, 1.f, 0.f};
  const vector3 mom{0.f, 0.f, 1.f};
  const detail::ray<test_algebra> r(pos, 0.f, mom, 0.f);

  // The bounding box
  mask<cuboid3D, test_algebra> c3{0u, x_min, y_min, z_min, x_max, y_max, z_max};
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{c3, 0u, envelope};

  ASSERT_TRUE(aabb.intersect(r));
}
