// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/geometry/shapes/annulus2D.hpp"

#include "detray/definitions/math.hpp"
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

constexpr scalar tol{1e-5f};

/// This tests the basic functionality of a stereo annulus
GTEST_TEST(detray_masks, annulus2D) {
  static_assert(concepts::shape<annulus2D, test_algebra>);
  static_assert(concepts::planar_shape<annulus2D, test_algebra>);

  constexpr scalar minR{7.2f * unit<scalar>::mm};
  constexpr scalar maxR{12.0f * unit<scalar>::mm};
  constexpr scalar minPhi{0.74195f};
  constexpr scalar maxPhi{1.33970f};
  point3 offset = {-2.f, 2.f, 0.f};

  // points in cartesian module frame
  point3 p2_in = {7.f, 7.f, 0.f};
  point3 p2_out1 = {5.f, 5.f, 0.f};
  point3 p2_out2 = {10.f, 3.f, 0.f};
  point3 p2_out3 = {10.f, 10.f, 0.f};
  point3 p2_out4 = {4.f, 10.f, 0.f};

  auto toStripFrame = [&offset](const point3 &xy) -> point3 {
    auto shifted = xy + offset;
    scalar r{vector::perp(shifted)};
    scalar phi{vector::phi(shifted)};
    return point3{r, phi, static_cast<scalar>(0.f)};
  };

  mask<annulus2D, test_algebra> ann2{0u,     minR, maxR,      minPhi,
                                     maxPhi, 0.f,  offset[0], offset[1]};

  ASSERT_NEAR(ann2[annulus2D::e_min_r], 7.2f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_max_r], 12.0f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_min_phi_rel], 0.74195f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_max_phi_rel], 1.33970f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_shift_x], -2.0f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_shift_y], 2.0f, tol);
  ASSERT_NEAR(ann2[annulus2D::e_average_phi], 0.f, tol);

  ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_in)));
  ASSERT_FALSE(ann2.is_inside(toStripFrame(p2_out1)));
  ASSERT_FALSE(ann2.is_inside(toStripFrame(p2_out2)));
  ASSERT_FALSE(ann2.is_inside(toStripFrame(p2_out3)));
  ASSERT_FALSE(ann2.is_inside(toStripFrame(p2_out4)));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1), 1.3f));
  ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4), 0.9f));

  ASSERT_NEAR(
      ann2.get_shape().min_dist_to_boundary(ann2.values(), toStripFrame(p2_in)),
      2.1005f, tol);
  ASSERT_NEAR(ann2.get_shape().min_dist_to_boundary(ann2.values(),
                                                    toStripFrame(p2_out1)),
              0.128932f, tol);
  ASSERT_NEAR(ann2.get_shape().min_dist_to_boundary(ann2.values(),
                                                    toStripFrame(p2_out2)),
              1.55969f, tol);
  ASSERT_NEAR(ann2.get_shape().min_dist_to_boundary(ann2.values(),
                                                    toStripFrame(p2_out3)),
              2.14214f, tol);
  ASSERT_NEAR(ann2.get_shape().min_dist_to_boundary(ann2.values(),
                                                    toStripFrame(p2_out4)),
              0.80214f, tol);

  // Check area: @TODO not implemented, yet
  scalar a = ann2.area();
  ASSERT_EQ(a, ann2.measure());

  // Check corner positions
  std::array<scalar, 8> c = ann2.get_shape().corners(ann2.values());
  for (unsigned int i{0u}; i < 8u; i += 2u) {
    // Transform to local cartesian beam system
    const scalar loc_x{c[i] * math::cos(c[i + 1]) -
                       ann2.values()[annulus2D::e_shift_x]};
    const scalar loc_y{c[i] * math::sin(c[i + 1]) -
                       ann2.values()[annulus2D::e_shift_y]};

    // Inner points
    if (i < 4u) {
      EXPECT_NEAR(std::hypot(loc_x, loc_y), minR, tol)
          << "point " << i << ": loc_x: " << loc_x << ", loc_y: " << loc_y
          << std::endl;
    }
    // Outer points
    else {
      EXPECT_NEAR(std::hypot(loc_x, loc_y), maxR, tol)
          << "point " << i << ": loc_x: " << loc_x << ", loc_y: " << loc_y
          << std::endl;
    }
  }

  // Check bounding box
  constexpr scalar envelope{0.01f};
  const auto loc_bounds = ann2.local_min_bounds(envelope);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_x], 3.8954f - envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_y], 2.39186f - envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_min_z], -envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_x], 10.50652f + envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_y], 10.89317f + envelope, tol);
  ASSERT_NEAR(loc_bounds[cuboid3D::e_max_z], envelope, tol);

  // TODO: Check centroid against visualization
}

/// This tests the inside/outside method of the mask
GTEST_TEST(detray_masks, annulus2D_ratio_test) {
  struct mask_check {
    bool operator()(const point3 &p, const mask<annulus2D, test_algebra> &ann,
                    const test::transform3 &trf, const scalar t) {
      return ann.is_inside(trf, p, t);
    }
  };

  constexpr scalar t{0.f};

  constexpr mask<annulus2D, test_algebra> ann{0u,       2.5f, 5.f,  -0.64299f,
                                              4.13173f, 1.f,  0.5f, 0.f};
  const test::transform3 trf{};

  constexpr scalar world{10.f * unit<scalar>::mm};
  const auto n_points{static_cast<std::size_t>(std::pow(500, 3))};

  // x- and y-coordinates yield a valid local position on the underlying plane
  std::vector<point3> points =
      test::generate_regular_points<cuboid3D>(n_points, {world});

  scalar ratio = test::ratio_test<mask_check>(points, ann, trf, t);

  ASSERT_FALSE(detail::is_invalid_value(ratio));
}
