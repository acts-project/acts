// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/utils/curvilinear_frame.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using transform3 = test::transform3;
using vector3 = typename transform3::vector3;
using point3 = test::point3;
using scalar = test::scalar;

GTEST_TEST(detray_utils, curvilinear_frame) {
  constexpr const scalar tolerance = 1e-5f;

  point3 pos = {4.f, 10.f, 2.f};
  scalar time = 0.1f;
  vector3 mom = {10.f, 20.f, 30.f};
  constexpr scalar charge = -1.f;

  free_track_parameters<test_algebra> free_params(pos, time, mom, charge);

  curvilinear_frame<test_algebra> cf(free_params);

  const auto bound_vec = cf.m_bound_vec;

  const scalar phi = vector::phi(mom);
  const scalar theta = vector::theta(mom);

  EXPECT_NEAR(getter::element(bound_vec.bound_local(), e_bound_loc0), 0.f,
              tolerance);
  EXPECT_NEAR(getter::element(bound_vec.bound_local(), e_bound_loc1), 0.f,
              tolerance);
  EXPECT_NEAR(bound_vec.phi(), phi, tolerance);
  EXPECT_NEAR(bound_vec.theta(), theta, tolerance);
  EXPECT_NEAR(bound_vec.qop(), charge / vector::norm(mom), tolerance);

  const vector3 unit_p = vector::normalize(mom);
  const vector3 trf_x = cf.m_trf.x();
  const vector3 trf_y = cf.m_trf.y();
  const vector3 trf_z = cf.m_trf.z();
  const vector3 trf_t = cf.m_trf.translation();

  // local x axis transform does not have global z components
  EXPECT_NEAR(trf_x[2], 0.f, tolerance);

  // local z axis of transform = track direction
  EXPECT_NEAR(unit_p[0], trf_z[0], tolerance);
  EXPECT_NEAR(unit_p[1], trf_z[1], tolerance);
  EXPECT_NEAR(unit_p[2], trf_z[2], tolerance);

  // translation of transform = track position
  EXPECT_NEAR(pos[0], trf_t[0], tolerance);
  EXPECT_NEAR(pos[1], trf_t[1], tolerance);
  EXPECT_NEAR(pos[2], trf_t[2], tolerance);

  // Top-left components of the jacobian
  const auto jac = cf.bound_to_free_jacobian();
  EXPECT_NEAR(getter::element(jac, 0u, 0u), trf_x[0u], tolerance);
  EXPECT_NEAR(getter::element(jac, 1u, 0u), trf_x[1u], tolerance);
  EXPECT_NEAR(getter::element(jac, 2u, 0u), trf_x[2u], tolerance);
  EXPECT_NEAR(getter::element(jac, 0u, 1u), trf_y[0u], tolerance);
  EXPECT_NEAR(getter::element(jac, 1u, 1u), trf_y[1u], tolerance);
  EXPECT_NEAR(getter::element(jac, 2u, 1u), trf_y[2u], tolerance);

  EXPECT_NEAR(getter::element(jac, e_free_dir0, e_bound_phi),
              -math::sin(phi) * math::sin(theta), tolerance);
  EXPECT_NEAR(getter::element(jac, e_free_dir0, e_bound_theta),
              math::cos(phi) * math::cos(theta), tolerance);
  EXPECT_NEAR(getter::element(jac, e_free_dir1, e_bound_phi),
              math::cos(phi) * math::sin(theta), tolerance);
  EXPECT_NEAR(getter::element(jac, e_free_dir1, e_bound_theta),
              math::sin(phi) * math::cos(theta), tolerance);
  EXPECT_NEAR(getter::element(jac, e_free_dir2, e_bound_phi), 0, tolerance);
  EXPECT_NEAR(getter::element(jac, e_free_dir2, e_bound_theta),
              -math::sin(theta), tolerance);
}
