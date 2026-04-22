// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/bounding_volume.hpp"

#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google test include(s).
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;

namespace {

// envelope around wrapped object (scalor in percent)
constexpr scalar envelope{0.01f};

// test tolerance
constexpr scalar tol{1e-5f};

}  // anonymous namespace

/// This tests the basic functionality cuboid axis aligned bounding box
GTEST_TEST(detray_utils, bounding_cuboid3D) {
  // cuboid
  constexpr scalar hx{1.f * unit<scalar>::mm};
  constexpr scalar hy{9.3f * unit<scalar>::mm};
  constexpr scalar hz{0.5f * unit<scalar>::mm};

  point3 p2_in = {0.5f, 8.0f, -0.4f};
  point3 p2_edge = {1.f, 9.3f, 0.5f};
  point3 p2_out = {1.5f, -9.f, 0.55f};

  mask<cuboid3D, test_algebra> c3{0u, -hx, -hy, -hz, hx, hy, hz};

  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{c3, 0u, envelope};

  // Id of this instance
  ASSERT_EQ(aabb.id(), 0u);

  // Test the bounds
  const auto bounds = aabb.bounds();
  ASSERT_NEAR(bounds[cuboid3D::e_min_x], -hx - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_min_y], -hy - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_min_z], -hz - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_x], hx + envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_y], hy + envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_z], hz + envelope, tol);

  ASSERT_TRUE(aabb.is_inside(p2_in));
  ASSERT_TRUE(aabb.is_inside(p2_edge));
  ASSERT_FALSE(aabb.is_inside(p2_out));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(aabb.is_inside(p2_out, 1.f));
}

/// This tests the basic functionality cylindrical axis aligned bounding box
GTEST_TEST(detray_utils, bounding_cylinder3D) {
  // cylinder
  constexpr scalar r{3.f * unit<scalar>::mm};
  constexpr scalar hz{4.f * unit<scalar>::mm};

  point3 p3_in = {r, static_cast<scalar>(0.f), static_cast<scalar>(-1.f)};
  point3 p3_edge = {static_cast<scalar>(0.f), r, hz};
  point3 p3_out = {r * constant<scalar>::inv_sqrt2,
                   r * constant<scalar>::inv_sqrt2, static_cast<scalar>(4.5f)};
  point3 p3_off = {1.f, 1.f, -9.f};

  axis_aligned_bounding_volume<cylinder3D, test_algebra> aabc{
      0u, 0.f, -constant<scalar>::pi, -hz, r, constant<scalar>::pi, hz};

  // Id of this instance
  ASSERT_EQ(aabc.id(), 0u);

  // Test the bounds
  const auto bounds = aabc.bounds();
  ASSERT_NEAR(bounds[cylinder3D::e_min_r], 0.f, tol);
  ASSERT_NEAR(bounds[cylinder3D::e_max_r], r, tol);
  ASSERT_NEAR(bounds[cylinder3D::e_min_phi], -constant<scalar>::pi, tol);
  ASSERT_NEAR(bounds[cylinder3D::e_max_phi], constant<scalar>::pi, tol);
  ASSERT_NEAR(bounds[cylinder3D::e_min_z], -hz, tol);
  ASSERT_NEAR(bounds[cylinder3D::e_max_z], hz, tol);

  ASSERT_TRUE(aabc.is_inside(p3_in));
  ASSERT_TRUE(aabc.is_inside(p3_edge));
  ASSERT_FALSE(aabc.is_inside(p3_out));
  ASSERT_FALSE(aabc.is_inside(p3_off));
  // Move outside point inside using a tolerance
  ASSERT_TRUE(aabc.is_inside(p3_out, 0.6f));
}

/// This tests the basic functionality of an aabb around a stereo annulus
GTEST_TEST(detray_utils, annulus2D_aabb) {
  constexpr scalar minR{7.2f * unit<scalar>::mm};
  constexpr scalar maxR{12.0f * unit<scalar>::mm};
  constexpr scalar minPhi{0.74195f};
  constexpr scalar maxPhi{1.33970f};
  typename transform3::point2 offset = {-2.f, 2.f};

  mask<annulus2D, test_algebra> ann2{0u,     minR, maxR,      minPhi,
                                     maxPhi, 0.f,  offset[0], offset[1]};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{ann2, 0u, envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  const auto glob_aabb = aabb.transform(trf);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_x], 2.39186f - envelope + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_y], -10.50652f - envelope + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_z], -envelope + t[2], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_x], 10.89317f + envelope + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_y], -3.8954f + envelope + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_z], envelope + t[2], tol);
}

/// This tests the basic functionality of an aabb around a cylinder
GTEST_TEST(detray_utils, cylinder2D_aabb) {
  constexpr scalar r{3.f * unit<scalar>::mm};
  constexpr scalar hz{4.f * unit<scalar>::mm};

  mask<cylinder2D, test_algebra> c{0u, r, -hz, hz};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{c, 0u, envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  const auto glob_aabb = aabb.transform(trf);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_x], -(r + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_y], -(r + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_z], -(hz + envelope) + t[2], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_x], (r + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_y], (r + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_z], (hz + envelope) + t[2], tol);
}

/// This tests the basic functionality of an aabb around a line
GTEST_TEST(detray_utils, line2D_aabb) {
  constexpr scalar cell_size{1.f * unit<scalar>::mm};
  constexpr scalar hz{50.f * unit<scalar>::mm};

  const mask<line_circular, test_algebra> ln_r{0u, cell_size, hz};
  const mask<line_square, test_algebra> ln_sq{0u, cell_size, hz};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb_r{ln_r, 0u,
                                                              envelope};
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb_sq{ln_sq, 0u,
                                                               envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  // Radial crossection
  const auto glob_aabb_r = aabb_r.transform(trf);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_min_x], -(cell_size + envelope) + t[0],
              tol);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_min_y], -(cell_size + envelope) + t[1],
              tol);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_min_z], -(hz + envelope) + t[2], tol);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_max_x], (cell_size + envelope) + t[0],
              tol);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_max_y], (cell_size + envelope) + t[1],
              tol);
  ASSERT_NEAR(glob_aabb_r[cuboid3D::e_max_z], (hz + envelope) + t[2], tol);

  // Square crossection
  const auto glob_aabb_sq = aabb_sq.transform(trf);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_min_x], -(cell_size + envelope) + t[0],
              tol);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_min_y], -(cell_size + envelope) + t[1],
              tol);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_min_z], -(hz + envelope) + t[2], tol);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_max_x], (cell_size + envelope) + t[0],
              tol);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_max_y], (cell_size + envelope) + t[1],
              tol);
  ASSERT_NEAR(glob_aabb_sq[cuboid3D::e_max_z], (hz + envelope) + t[2], tol);
}

/// This tests the basic functionality of an aabb around a rectangle
GTEST_TEST(detray_utils, rectangle2D_aabb) {
  constexpr scalar hx{1.f * unit<scalar>::mm};
  constexpr scalar hy{9.3f * unit<scalar>::mm};

  mask<rectangle2D, test_algebra> r2{0u, hx, hy};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{r2, 0u, envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  const auto glob_aabb = aabb.transform(trf);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_x], -(hy + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_y], -(hx + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_z], -(envelope) + t[2], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_x], (hy + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_y], (hx + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_z], (envelope) + t[2], tol);
}

/// This tests the basic functionality of an aabb around a rectangle
GTEST_TEST(detray_utils, ring2D_aabb) {
  constexpr scalar inner_r{0.f * unit<scalar>::mm};
  constexpr scalar outer_r{3.5f * unit<scalar>::mm};

  mask<ring2D, test_algebra> r2{0u, inner_r, outer_r};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{r2, 0u, envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  const auto glob_aabb = aabb.transform(trf);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_x], -(outer_r + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_y], -(outer_r + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_z], -(envelope) + t[2], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_x], (outer_r + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_y], (outer_r + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_z], (envelope) + t[2], tol);
}

/// This tests the basic functionality of an aabb around a rectangle
GTEST_TEST(detray_utils, trapezoid2D_aabb) {
  constexpr scalar hx_miny{1.f * unit<scalar>::mm};
  constexpr scalar hx_maxy{3.f * unit<scalar>::mm};
  constexpr scalar hy{2.f * unit<scalar>::mm};
  constexpr scalar divisor{1.f / (2.f * hy)};

  mask<trapezoid2D, test_algebra> t2{0u, hx_miny, hx_maxy, hy, divisor};

  // Construct local aabb around mask
  axis_aligned_bounding_volume<cuboid3D, test_algebra> aabb{t2, 0u, envelope};

  // rotate around z-axis by 90deg and then translate by 1mm in each direction
  const vector3 new_x{0.f, -1.f, 0.f};
  const vector3 new_z{0.f, 0.f, 1.f};
  const vector3 t{1.f, 2.f, 3.f};
  const transform3 trf{t, new_z, new_x};

  const auto glob_aabb = aabb.transform(trf);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_x], -(hy + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_y], -(hx_maxy + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_min_z], -envelope + t[2], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_x], (hy + envelope) + t[0], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_y], (hx_maxy + envelope) + t[1], tol);
  ASSERT_NEAR(glob_aabb[cuboid3D::e_max_z], envelope + t[2], tol);
}

/// This tests wrapping a collection of cuboid bounding volumes
GTEST_TEST(detray_utils, wrap_bounding_cuboid3D) {
  using box_t = axis_aligned_bounding_volume<cuboid3D, test_algebra>;

  box_t b1{0u, -1.f, 0.f, -10.f, -0.5f, 1.f, 0.f};
  box_t b2{0u, -2.f, 3.f, 2.f, 2.f, 4.5f, 3.5f};
  box_t b3{0u, -1.5f, -0.1f, -2.f, 0.f, 0.1f, 2.f};
  box_t b4{0u, 0.f, 0.f, 0.f, 2.f, 2.f, 10.f};

  std::vector<box_t> boxes = {b1, b2, b3, b4};

  box_t aabb{boxes, 0u, envelope};

  // Id of this instance
  ASSERT_EQ(aabb.id(), 0u);

  // Test the bounds
  const auto bounds = aabb.bounds();
  ASSERT_NEAR(bounds[cuboid3D::e_min_x], -2.f - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_min_y], -0.1f - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_min_z], -10.f - envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_x], 2.f + envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_y], 4.5f + envelope, tol);
  ASSERT_NEAR(bounds[cuboid3D::e_max_z], 10.f + envelope, tol);
}
