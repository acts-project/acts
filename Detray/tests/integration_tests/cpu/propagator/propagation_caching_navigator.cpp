// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tracks/trajectories.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;

constexpr scalar tol{1e-3f};
constexpr scalar path_limit{5.f * unit<scalar>::cm};
constexpr std::size_t cache_size{navigation::default_cache_size};

namespace detray::actor {

/// Compare helical track positions for stepper
struct helix_inspector : public base_actor {
  /// Keeps the state of a helix gun to calculate track positions
  struct state {
    scalar path_from_surface{0.f};
  };

  /// Check that the stepper remains on the right helical track for its pos.
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(
      state& inspector_state, const propagator_state_t& prop_state,
      parameter_transporter_result<test_algebra>& res) const {
    const auto& navigation = prop_state.navigation();
    const auto& stepping = prop_state.stepping();
    const auto& bound_params = res.destination_params();

    typename propagator_state_t::detector_type::geometry_context ctx{};

    // Nothing has happened yet (first call of actor chain)
    if (stepping.path_length() < tol ||
        inspector_state.path_from_surface < tol) {
      return;
    }

    if (bound_params.surface_link().is_invalid()) {
      return;
    }

    // Surface
    const auto sf =
        tracking_surface{navigation.detector(), bound_params.surface_link()};

    const free_track_parameters<test_algebra> free_params =
        sf.bound_to_free_vector(ctx, bound_params);

    const auto last_pos = free_params.pos();

    const auto bvec =
        stepping.field().at(last_pos[0], last_pos[1], last_pos[2]);
    const vector3 b{bvec[0], bvec[1], bvec[2]};

    detail::helix<test_algebra> hlx(free_params, b);

    const auto true_pos = hlx(inspector_state.path_from_surface);

    const point3 relative_error{1.f / inspector_state.path_from_surface *
                                (stepping().pos() - true_pos)};

    ASSERT_NEAR(vector::norm(relative_error), 0.f, tol);

    auto true_J = hlx.jacobian(inspector_state.path_from_surface);

    for (unsigned int i = 0u; i < e_free_size; i++) {
      for (unsigned int j = 0u; j < e_free_size; j++) {
        ASSERT_NEAR(getter::element(stepping.transport_jacobian(), i, j),
                    getter::element(true_J, i, j),
                    inspector_state.path_from_surface * tol * 10.f);
      }
    }
    // The propagation does not start on a surface, skip the initial path
    if (!bound_params.surface_link().is_invalid()) {
      inspector_state.path_from_surface += stepping.step_size();
    }
    // Reset path from surface
    if (navigation.is_on_sensitive() || navigation.encountered_sf_material()) {
      inspector_state.path_from_surface = 0.f;
    }
  }
};

}  // namespace detray::actor

/// Test basic functionality of the propagator using a straight line stepper
GTEST_TEST(detray_propagator, propagator_line_stepper) {
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(false);
  const auto [d, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  using navigator_t =
      caching_navigator<decltype(d), cache_size, navigation::print_inspector>;
  using stepper_t = line_stepper<test_algebra>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

  const point3 pos{0.f, 0.f, 0.f};
  const vector3 mom{1.f, 1.f, 0.f};
  free_track_parameters<test_algebra> track(pos, 0.f, mom, -1.f);

  propagation::config prop_cfg{};
  propagator_t p{prop_cfg};

  propagator_t::state state(track, d, prop_cfg.context);

  EXPECT_TRUE(p.propagate(state))
      << state.navigation().inspector().to_string() << std::endl;
}

/// Fixture for Runge-Kutta Propagation
class PropagatorWithRkStepper : public ::testing::TestWithParam<
                                    std::tuple<scalar, scalar, test::vector3>> {
 public:
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  /// Set the test environment up
  virtual void SetUp() {
    overstep_tol = std::get<0>(GetParam());
    step_constr = std::get<1>(GetParam());

    trk_gen_cfg.phi_steps(50u).theta_steps(50u);
    trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);
  }

  /// Clean up
  virtual void TearDown() { /* Do nothing */ }

 protected:
  /// Detector configuration
  vecmem::host_memory_resource host_mr;

  /// Toy detector configuration
  toy_det_config<scalar> toy_cfg =
      toy_det_config<scalar>{}.n_brl_layers(4u).n_edc_layers(7u);

  /// Track generator configuration
  generator_t::configuration trk_gen_cfg{};

  /// Stepper configuration
  scalar overstep_tol;
  scalar step_constr;
};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_const_bfield) {
  // Constant magnetic field type
  using bfield_t = bfield::const_field_t<scalar>;

  // Toy detector
  using detector_t = detector<test::toy_metadata>;

  // Runge-Kutta propagation
  using navigator_t =
      caching_navigator<detector_t, cache_size, navigation::print_inspector>;
  using track_t = free_track_parameters<test_algebra>;
  using constraints_t = constrained_step<scalar>;
  using policy_t = stepper_rk_policy<scalar>;
  using stepper_t =
      rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
  // Include helix actor to check track position/covariance
  using actor_chain_t = actor_chain<
      actor::pathlimit_aborter<scalar>,
      actor::parameter_updater<
          test_algebra, actor::pointwise_material_interactor<test_algebra>,
          actor::helix_inspector>>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Build detector
  toy_cfg.use_material_maps(false);
  toy_cfg.mapped_material(detray::vacuum<scalar>());
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  const bfield_t bfield = create_const_field<scalar>(std::get<2>(GetParam()));

  // Propagator is built from the stepper and navigator
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(overstep_tol);
  cfg.navigation.search_window = {3u, 3u};
  cfg.navigation.estimate_scattering_noise = false;
  propagator_t p{cfg};

  // Iterate through uniformly distributed momentum directions
  for (auto track : generator_t{trk_gen_cfg}) {
    assert(track.qop() != 0.f);

    // Generate second track state used for propagation with pathlimit
    track_t lim_track(track);

    // Build actor states: the helix inspector can be shared
    auto actor_states = actor_chain_t::make_default_actor_states();
    auto actor_states_lim = actor_chain_t::make_default_actor_states();

    // Make sure the lim state is being terminated
    auto& pathlimit_aborter_state =
        detail::get<actor::pathlimit_aborter<scalar>::state>(actor_states_lim);
    pathlimit_aborter_state.set_path_limit(path_limit);

    // Init propagator states
    typename propagator_t::stepper_type::magnetic_field_type bfield_view(
        bfield);
    propagator_t::state state(track, bfield_view, det);
    propagator_t::state sync_state(track, bfield_view, det);
    propagator_t::state lim_state(lim_track, bfield_view, det);

    state.debug(true);
    sync_state.debug(true);
    lim_state.debug(true);

    // Set step constraints
    state.stepping().template set_constraint<step::constraint::e_accuracy>(
        step_constr);
    sync_state.stepping().template set_constraint<step::constraint::e_accuracy>(
        step_constr);
    lim_state.stepping().template set_constraint<step::constraint::e_accuracy>(
        step_constr);

    // No multiple scattering simulated in this test
    using updater_state_t = actor::parameter_updater_state<test_algebra>;
    detail::get<updater_state_t>(actor_states)
        .noise_estimation_cfg()
        .estimate_scattering_noise = false;
    detail::get<updater_state_t>(actor_states_lim)
        .noise_estimation_cfg()
        .estimate_scattering_noise = false;

    // Propagate the entire detector
    ASSERT_TRUE(
        p.propagate(state, actor_chain_t::setup_actor_states(actor_states)))
        << state.navigation().inspector().to_string() << std::endl;

    // Propagate with path limit
    ASSERT_FALSE(p.propagate(
        lim_state, actor_chain_t::setup_actor_states(actor_states_lim)))
        << lim_state.navigation().inspector().to_string() << std::endl;

    ASSERT_GE(std::abs(path_limit), lim_state.stepping().abs_path_length())
        << "Absolute path length: " << lim_state.stepping().abs_path_length()
        << ", path limit: " << path_limit << std::endl;
    //<< state.navigation().inspector().to_string() << std::endl;
  }
}

/// Test propagation in an inhomogeneous magnetic field using a Runge-Kutta
/// stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_inhom_bfield) {
  // Magnetic field map using nearest neighbor interpolation
  using bfield_t = bfield::inhom_field_t<scalar>;

  // Toy detector
  using detector_t = detector<test::toy_metadata>;

  // Runge-Kutta propagation
  using navigator_t =
      caching_navigator<detector_t, cache_size, navigation::print_inspector>;
  using track_t = free_track_parameters<test_algebra>;
  using constraints_t = constrained_step<scalar>;
  using policy_t = stepper_rk_policy<scalar>;
  using stepper_t =
      rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
  // Include helix actor to check track position/covariance
  using actor_chain_t = actor_chain<
      actor::pathlimit_aborter<scalar>,
      actor::parameter_updater<
          test_algebra, actor::pointwise_material_interactor<test_algebra>>>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Build detector and magnetic field
  toy_cfg.use_material_maps(false);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);
  const bfield_t bfield = create_inhom_field<scalar>();

  // Propagator is built from the stepper and navigator
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(overstep_tol);
  cfg.navigation.search_window = {3u, 3u};
  cfg.navigation.estimate_scattering_noise = false;
  propagator_t p{cfg};

  // Iterate through uniformly distributed momentum directions
  for (auto track : generator_t{trk_gen_cfg}) {
    // Generate second track state used for propagation with pathlimit
    track_t lim_track(track);

    // Build actor states: the helix inspector can be shared
    actor::pathlimit_aborter<scalar>::state unlimted_aborter_state{};
    actor::pathlimit_aborter<scalar>::state pathlimit_aborter_state{path_limit};
    actor::parameter_updater_state<test_algebra> updater_state{cfg};
    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    // Create actor states tuples
    auto actor_states =
        detray::tie(unlimted_aborter_state, updater_state, interactor_state);
    auto lim_actor_states =
        detray::tie(pathlimit_aborter_state, updater_state, interactor_state);

    // Init propagator states
    typename propagator_t::stepper_type::magnetic_field_type bfield_view(
        bfield);
    propagator_t::state state(track, bfield_view, det);
    propagator_t::state lim_state(lim_track, bfield_view, det);

    // Set step constraints
    state.stepping().template set_constraint<step::constraint::e_accuracy>(
        step_constr);
    lim_state.stepping().template set_constraint<step::constraint::e_accuracy>(
        step_constr);

    // Propagate the entire detector
    state.debug(true);
    ASSERT_TRUE(p.propagate(state, actor_states))
        << state.navigation().inspector().to_string() << std::endl;

    // Propagate with path limit
    ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, tol);
    lim_state.debug(true);
    ASSERT_FALSE(p.propagate(lim_state, lim_actor_states))
        << lim_state.navigation().inspector().to_string() << std::endl;

    ASSERT_TRUE(lim_state.stepping().path_length() < std::abs(path_limit) + tol)
        << "path length: " << lim_state.stepping().path_length()
        << ", path limit: " << path_limit << std::endl;
    //<< state.navigation().inspector().to_string() << std::endl;
  }
}

// No step size constraint
INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation1, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-100.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              2.f * unit<scalar>::T})));

// Add some restrictions for more frequent navigation updates in the cases of
// non-z-aligned B-fields
INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation2, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-400.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation3, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-400.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation4, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-600.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));
