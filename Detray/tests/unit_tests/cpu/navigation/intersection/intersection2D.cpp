// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/invalid_values.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

using scalar = test::scalar;

namespace {

/// Define mask types
enum class mask_id : unsigned int {
  e_unmasked = 0u,
};

/// Define material types
enum class material_id : unsigned int {
  e_material_slab = 0u,
};

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

}  // namespace

using test_algebra = test::algebra;
using scalar_t = dscalar<test_algebra>;
using point2 = test::point2;
using vector3 = test::vector3;
using point3 = test::point3;

using mask_link_t = dtyped_index<mask_id, dindex>;
using material_link_t = dtyped_index<material_id, dindex>;
using surface_t =
    surface_descriptor<mask_link_t, material_link_t, test_algebra>;

// This tests the construction of a intresection
GTEST_TEST(detray_intersection, intersection2D) {
  using enum detray::intersection::status;

  using intersection_t =
      intersection2D<surface_t, test_algebra, intersection::contains_pos>;
  using nominal_inters_t =
      intersection2D<surface_t, test_algebra, !intersection::contains_pos>;

  // Check memory layout of intersection struct (uncomment for debugging)
  /*static_assert(offsetof(nominal_inters_t, m_surface) == 0);
  static_assert(offsetof(nominal_inters_t, m_ip) == 16);
  // Depends on floating point precision of 'path' member variable
  static_assert((offsetof(nominal_inters_t, m_volume_link) == 20) ||
                (offsetof(nominal_inters_t, m_volume_link) == 24));
  static_assert((offsetof(nominal_inters_t, m_status) == 22) ||
                (offsetof(nominal_inters_t, m_status) == 26));
  static_assert((offsetof(nominal_inters_t, m_direction) == 23) ||
                (offsetof(nominal_inters_t, m_direction) == 27));*/

  // 24 bytes for single precision, 32 bytes for double
  static_assert((sizeof(nominal_inters_t) == 24) ||
                (sizeof(nominal_inters_t) == 32));

  const surface_t sf{};
  const point3 test_pt{0.2f, 0.4f, 0.f};

  intersection_t i0{sf, 2.f, test_pt, 1u, e_outside, true};
  intersection_t i1{sf, 1.7f, test_pt, 0u, e_inside, false};

  intersection_t invalid{};
  ASSERT_FALSE(invalid.is_inside());

  dvector<intersection_t> intersections = {invalid, i0, i1};
  std::ranges::sort(intersections);

  ASSERT_NEAR(intersections[0].path(), 1.7f, tol);
  ASSERT_NEAR(intersections[1].path(), 2.f, tol);
  ASSERT_TRUE(detail::is_invalid_value(intersections[2].path()));
}
