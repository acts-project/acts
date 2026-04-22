// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/geometry/shapes/line.hpp"

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

namespace {

constexpr scalar tol{1e-4f};

// 50 mm wire with 1 mm radial cell size
constexpr scalar cell_size{1.f * unit<scalar>::mm};
constexpr scalar hz{50.f * unit<scalar>::mm};

}  // anonymous namespace

/// This tests the basic functionality of a line with a radial cross section
GTEST_TEST(detray_masks, line_circular) {
  static_assert(concepts::shape<line_circular, test_algebra>);
  static_assert(concepts::line_shape<line_circular, test_algebra>);

  const point3 ln_in{0.09f, 0.5f, 0.f};
  const point3 ln_edge{1.f, 50.f, 0.f};
  const point3 ln_out1{1.2f, 0.f, 0.f};
  const point3 ln_out2{0.09f, -51.f, 0.f};

  const mask<line_circular, test_algebra> ln{0u, cell_size, hz};

  ASSERT_NEAR(ln[line_circular::e_cross_section], 1.f * unit<scalar>::mm, tol);
  ASSERT_NEAR(ln[line_circular::e_half_z], 50.f * unit<scalar>::mm, tol);

  ASSERT_TRUE(ln.is_inside(ln_in));
  ASSERT_TRUE(ln.is_inside(ln_edge));
  ASSERT_FALSE(ln.is_inside(ln_out1));
  ASSERT_FALSE(ln.is_inside(ln_out2));

  // Check area and measure
  EXPECT_NEAR(ln.area(), 628.318531f * unit<scalar>::mm2, tol);
  EXPECT_NEAR(ln.measure(), 314.159265359f * unit<scalar>::mm3, tol);

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = ln.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -(hz + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], (hz + envelope), tol);

  const auto centroid = ln.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, line_circular_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p,
                    const mask<line_circular, test_algebra> &st,
                    const test::transform3 &trf, const test::vector3 & /*dir*/,
                    const scalar t) {
      return st.is_inside(trf, p, t);
    }
  };

  constexpr mask<line_circular, test_algebra> st{0u, 3.f, 5.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  // Track direction not parallel to the line, so that we always
  // have a valid local position (normally ensure by the intersector)
  test::vector3 dir{1.f, 0.f, 0.f};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, st, trf, dir, t);

  const scalar volume{st.measure()};
  const scalar world{size * size * size};

  ASSERT_NEAR(ratio, volume / world, 0.001f);
}

/// This tests the basic functionality of a line with a square cross section
GTEST_TEST(detray_masks, line_square) {
  static_assert(concepts::shape<line_square, test_algebra>);
  static_assert(concepts::line_shape<line_circular, test_algebra>);

  const point3 ln_in{static_cast<scalar>(1.1f), static_cast<scalar>(0.9f),
                     constant<scalar>::pi_4};
  const point3 ln_edge{1.f, 1.f, 0.f};
  const point3 ln_out1{1.1f, 0.f, 0.f};
  const point3 ln_out2{0.09f, -51.f, 0.f};

  // 50 mm wire with 1 mm square cell sizes
  const mask<line_square, test_algebra> ln{0u, cell_size, hz};

  ASSERT_NEAR(ln[line_circular::e_cross_section], 1.f * unit<scalar>::mm, tol);
  ASSERT_NEAR(ln[line_circular::e_half_z], 50.f * unit<scalar>::mm, tol);

  ASSERT_TRUE(ln.is_inside(ln_in));
  ASSERT_TRUE(ln.is_inside(ln_edge, 1e-5f));
  ASSERT_FALSE(ln.is_inside(ln_edge, -1e-5f));
  ASSERT_FALSE(ln.is_inside(ln_out1));
  ASSERT_FALSE(ln.is_inside(ln_out2));

  // Check area and measure
  EXPECT_NEAR(ln.area(), 800.f * unit<scalar>::mm2, tol);
  EXPECT_NEAR(ln.measure(), 400.f * unit<scalar>::mm3, tol);

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = ln.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], -(cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], -(cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -(hz + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], (cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], (cell_size + envelope), tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], (hz + envelope), tol);

  const auto centroid = ln.centroid();
  ASSERT_NEAR(centroid[0], 0.f, tol);
  ASSERT_NEAR(centroid[1], 0.f, tol);
  ASSERT_NEAR(centroid[2], 0.f, tol);
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, line_square_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<line_square, test_algebra> &dcl,
                    const test::transform3 &trf, const test::vector3 & /*dir*/,
                    const scalar t) {
      return dcl.is_inside(trf, p, t);
    }
  };

  constexpr mask<line_square, test_algebra> dcl{0u, 3.f, 5.f};

  constexpr scalar t{0.f};
  const test::transform3 trf{};
  // Track direction not parallel to the line, so that we always
  // have a valid local position (normally ensure by the intersector)
  test::vector3 dir{1.f, 0.f, 0.f};
  constexpr scalar size{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {size});

  scalar ratio = test::ratio_test<mask_check>(points, dcl, trf, dir, t);

  const scalar volume{dcl.measure()};
  const scalar world{size * size * size};

  ASSERT_NEAR(ratio, volume / world, 0.005f);
}
