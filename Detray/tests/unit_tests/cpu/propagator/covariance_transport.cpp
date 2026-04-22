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

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/shapes/unbounded.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/concepts.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/axis_rotation.hpp"

// Detray test include(s)
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

// Algebra types
using test_algebra = test::algebra;
using scalar = test::scalar;
using point2 = test::point2;
using vector3 = test::vector3;

// Mask types to be tested
// @TODO: Remove unbounded tag
using annulus_type =
    detray::mask<detray::unbounded<detray::annulus2D>, test_algebra>;
using rectangle_type = detray::mask<detray::rectangle2D, test_algebra>;
using trapezoid_type = detray::mask<detray::trapezoid2D, test_algebra>;
using ring_type = detray::mask<detray::ring2D, test_algebra>;
using cylinder_type = detray::mask<detray::cylinder2D, test_algebra>;
using straw_tube_type = detray::mask<detray::line_circular, test_algebra>;
using drift_cell_type = detray::mask<detray::line_square, test_algebra>;

constexpr scalar tol{1e-6f};
constexpr std::size_t cache_size{navigation::default_cache_size};

GTEST_TEST(detray_propagator, covariance_transport) {
  static_assert(
      detray::concepts::actor<actor::parameter_transporter<test_algebra>>);
  static_assert(detray::concepts::actor<actor::parameter_setter<test_algebra>>);

  vecmem::host_memory_resource host_mr;

  // Build in x-direction from given module positions
  detail::ray<test_algebra> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
  std::vector<scalar> positions = {0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f};

  tel_det_config<test_algebra, rectangle2D> tel_cfg{200.f * unit<scalar>::mm,
                                                    200.f * unit<scalar>::mm};
  tel_cfg.positions(positions).pilot_track(traj);

  // Build telescope detector with unbounded planes
  const auto [det, names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  using navigator_t =
      caching_navigator<decltype(det), cache_size, navigation::print_inspector>;
  using cline_stepper_t = line_stepper<test_algebra>;
  using actor_chain_t = actor_chain<actor::parameter_updater<test_algebra>>;
  using propagator_t = propagator<cline_stepper_t, navigator_t, actor_chain_t>;

  // Bound vector
  bound_parameters_vector<test_algebra> bound_vector{};
  bound_vector.set_theta(constant<scalar>::pi_4);
  bound_vector.set_qop(-0.1f);

  // Bound covariance
  auto bound_cov = matrix::identity<
      typename bound_track_parameters<test_algebra>::covariance_type>();

  // Note: Set angle error as ZERO, to constrain the loc0 and loc1 divergence
  getter::element(bound_cov, e_bound_phi, e_bound_phi) = 0.f;
  getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;

  // Bound track parameter
  const bound_track_parameters<test_algebra> bound_param0(
      det.surface(0u).identifier(), bound_vector, bound_cov);

  propagation::config prop_cfg{};
  prop_cfg.navigation.intersection.overstep_tolerance =
      -100.f * unit<float>::um;
  propagator_t p{prop_cfg};
  propagator_t::state propagation(bound_param0, det, prop_cfg.context);

  // Run propagator
  actor::parameter_updater_state<test_algebra> updater_state{prop_cfg,
                                                             bound_param0};
  updater_state.always_update();

  EXPECT_TRUE(p.propagate(propagation, detray::tie(updater_state)))
      << propagation.navigation().inspector().to_string();

  // Bound state after one turn propagation
  const auto& bound_param1 = updater_state.bound_params();

  // Check if the track reaches the final surface
  EXPECT_EQ(bound_param0.surface_link().volume(), 0u);
  EXPECT_EQ(bound_param0.surface_link().index(), 0u);
  EXPECT_EQ(bound_param1.surface_link().volume(), 0u);
  EXPECT_EQ(bound_param1.surface_link().id(), surface_id::e_sensitive);
  EXPECT_EQ(bound_param1.surface_link().index(), 6u);

  const auto bound_cov0 = bound_param0.covariance();
  const auto bound_cov1 = bound_param1.covariance();

  // Check covaraince
  for (unsigned int i = 0u; i < e_bound_size; i++) {
    for (unsigned int j = 0u; j < e_bound_size; j++) {
      EXPECT_NEAR(getter::element(bound_cov0, i, j),
                  getter::element(bound_cov1, i, j), tol);
    }
  }
}
