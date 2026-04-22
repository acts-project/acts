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

// detray include(s)
#include "detray/propagator/line_stepper.hpp"

#include "detray/definitions/units.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// google-test include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;
using point3 = test::point3;
using transform3 = test::transform3;

namespace {

constexpr scalar tol{1e-3f};
constexpr scalar step_size{1.f * unit<scalar>::mm};
constexpr stepping::config step_cfg{};

}  // namespace

// This tests the base functionality of the line stepper
GTEST_TEST(detray_propagator, line_stepper) {
  using namespace step;

  // Line stepper with and without constrained stepping
  using line_stepper_t = line_stepper<test_algebra>;
  using cline_stepper_t = line_stepper<test_algebra, constrained_step<scalar>>;

  point3 pos{0.f, 0.f, 0.f};
  vector3 mom{1.f, 1.f, 0.f};
  free_track_parameters<test_algebra> track(pos, 0.f, mom, -1.f);
  free_track_parameters<test_algebra> c_track(pos, 0.f, mom, -1.f);

  line_stepper_t l_stepper;
  cline_stepper_t cl_stepper;

  line_stepper_t::state l_state{track};
  cline_stepper_t::state cl_state{c_track};

  // Test the setting of step constraints
  cl_state.template set_constraint<constraint::e_accuracy>(10.f *
                                                           unit<scalar>::mm);
  cl_state.template set_constraint<constraint::e_actor>(2.f * unit<scalar>::mm);
  cl_state.template set_constraint<constraint::e_aborter>(5.f *
                                                          unit<scalar>::mm);
  cl_state.template set_constraint<constraint::e_user>(0.5f * unit<scalar>::mm);
  ASSERT_NEAR(cl_state.constraints().template size<>(), 0.5f * unit<scalar>::mm,
              tol);

  // Release all except 'actor', then set 'user' again
  cl_state.template release_step<constraint::e_accuracy>();
  cl_state.template release_step<constraint::e_aborter>();
  cl_state.template release_step<constraint::e_user>();
  ASSERT_NEAR(cl_state.constraints().template size<>(), 2.f * unit<scalar>::mm,
              tol);

  cl_state.template set_constraint<constraint::e_user>(0.5f * unit<scalar>::mm);
  ASSERT_NEAR(cl_state.constraints().template size<>(), 0.5f * unit<scalar>::mm,
              tol);

  // Run a few steps
  ASSERT_TRUE(l_stepper.step(step_size, l_state, step_cfg));
  // Step constraint to half step size
  ASSERT_TRUE(cl_stepper.step(step_size, cl_state, step_cfg));
  ASSERT_TRUE(cl_stepper.step(step_size, cl_state, step_cfg));

  track = l_state();
  ASSERT_NEAR(track.pos()[0], constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(track.pos()[1], constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(track.pos()[2], 0.f, tol);

  c_track = cl_state();
  ASSERT_NEAR(c_track.pos()[0], constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(c_track.pos()[1], constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(c_track.pos()[2], 0.f, tol);

  ASSERT_TRUE(l_stepper.step(step_size, l_state, step_cfg));

  track = l_state();
  ASSERT_NEAR(track.pos()[0], constant<scalar>::sqrt2, tol);
  ASSERT_NEAR(track.pos()[1], constant<scalar>::sqrt2, tol);
  ASSERT_NEAR(track.pos()[2], 0.f, tol);
}
