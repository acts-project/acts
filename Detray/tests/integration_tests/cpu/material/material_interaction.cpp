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
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/unbounded.hpp"
#include "detray/material/interaction.hpp"
#include "detray/material/material.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/test/common/bfield.hpp"

// Detray test include(s)
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/utils/random_scatterer.hpp"
#include "detray/test/utils/statistics.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using covariance_t =
    typename bound_track_parameters<test_algebra>::covariance_type;
using interactor_t = actor::pointwise_material_interactor<test_algebra>;

static_assert(detray::concepts::actor<interactor_t>);

// Test is done for muon
namespace {
pdg_particle ptc = muon<scalar>();
}

// Material interaction test with telescope Geometry
GTEST_TEST(detray_material, telescope_geometry_energy_loss) {
  vecmem::host_memory_resource host_mr;

  // Build in x-direction from given module positions
  detail::ray<test_algebra> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
  std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                   300.f, 350.f, 400.f, 450.f, 500.f};

  const auto mat = silicon_tml<scalar>();
  constexpr scalar thickness{0.17f * unit<scalar>::cm};

  tel_det_config<test_algebra, rectangle2D> tel_cfg{20.f * unit<scalar>::mm,
                                                    20.f * unit<scalar>::mm};
  tel_cfg.positions(positions)
      .pilot_track(traj)
      .module_material(mat)
      .mat_thickness(thickness);

  const auto [det, names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  using navigator_t = caching_navigator<decltype(det)>;
  using stepper_t = line_stepper<test_algebra>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar>;
  using actor_chain_t =
      actor_chain<pathlimit_aborter_t,
                  actor::parameter_updater<test_algebra, interactor_t>>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Propagator is built from the stepper and navigator
  propagation::config prop_cfg{};
  prop_cfg.navigation.intersection.overstep_tolerance =
      -100.f * unit<float>::um;
  propagator_t p{prop_cfg};

  constexpr scalar iniP{10.f * unit<scalar>::GeV};

  // Bound vector
  bound_parameters_vector<test_algebra> bound_vector{};
  bound_vector.set_theta(constant<scalar>::pi_2);
  bound_vector.set_qop(ptc.charge() / iniP);

  auto bound_cov = matrix::identity<covariance_t>();
  getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 0.f;

  // bound track parameter at first physical plane
  const bound_track_parameters<test_algebra> bound_param(
      det.surface(0u).identifier(), bound_vector, bound_cov);

  pathlimit_aborter_t::state aborter_state{};
  actor::parameter_updater_state<test_algebra> updater_state{prop_cfg,
                                                             bound_param};
  interactor_t::state interactor_state{};

  // Create actor states tuples
  auto actor_states =
      detray::tie(aborter_state, updater_state, interactor_state);

  propagator_t::state state(bound_param, det);
  state.debug(true);

  // Propagate the entire detector
  ASSERT_TRUE(p.propagate(state, actor_states));

  // new momentum
  const scalar newP{updater_state.bound_params().p(ptc.charge())};

  // mass
  const auto mass = ptc.mass();

  // new energy
  const scalar newE{std::hypot(newP, mass)};

  // Initial energy
  const scalar iniE{std::hypot(iniP, mass)};

  // New qop variance
  const scalar new_var_qop{
      getter::element(updater_state.bound_params().covariance(), e_bound_qoverp,
                      e_bound_qoverp)};

  // Interaction object
  interaction<scalar> I;

  // Zero incidence angle
  const scalar cos_inc_ang{1.f};

  // Same material used for default telescope detector
  material_slab<scalar> slab(mat, thickness);

  // Path segment in the material
  const scalar path_segment{slab.path_segment(cos_inc_ang)};

  // Expected Bethe Stopping power for telescope geometry is estimated
  // as (number of planes * energy loss per plane assuming 1 GeV muon).
  // It is not perfectly precise as the track loses its energy during
  // propagation. However, since the energy loss << the track momentum,
  // the assumption is not very bad
  const scalar dE{
      I.compute_energy_loss_bethe_bloch(path_segment, slab.get_material(), ptc,
                                        {ptc, ptc.charge() / iniP}) *
      static_cast<scalar>(positions.size())};

  // Check if the new energy after propagation is close enough to the
  // expected value
  EXPECT_NEAR(newE, iniE - dE, 1e-5f) << "dE = " << dE;

  const scalar sigma_qop{I.compute_energy_loss_landau_sigma_QOverP(
      path_segment, slab.get_material(), ptc, {ptc, ptc.charge() / iniP})};

  const scalar dvar_qop{sigma_qop * sigma_qop *
                        static_cast<scalar>(positions.size() - 1u)};

  EXPECT_NEAR(new_var_qop, dvar_qop, 1e-10f);

  // @todo: Validate the backward direction case as well?
}

// Material interaction test with telescope Geometry
GTEST_TEST(detray_material, telescope_geometry_scattering_angle) {
  vecmem::host_memory_resource host_mr;

  // Build in x-direction from given module positions
  detail::ray<test_algebra> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
  std::vector<scalar> positions = {0.f};

  // To make sure that someone won't put more planes than one by accident
  EXPECT_EQ(positions.size(), 1u);

  // Material
  const auto mat = silicon_tml<scalar>();
  const scalar thickness = 100.f * unit<scalar>::cm;

  // Create telescope geometry
  tel_det_config<test_algebra, rectangle2D> tel_cfg{2000.f * unit<scalar>::mm,
                                                    2000.f * unit<scalar>::mm};
  tel_cfg.positions(positions)
      .pilot_track(traj)
      .module_material(mat)
      .mat_thickness(thickness);

  const auto [det, names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  using navigator_t = caching_navigator<decltype(det)>;
  using stepper_t = line_stepper<test_algebra>;
  using simulator_t = actor::random_scatterer<test_algebra>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar>;
  using actor_chain_t =
      actor_chain<pathlimit_aborter_t,
                  actor::parameter_updater<test_algebra, simulator_t>>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Propagator is built from the stepper and navigator
  propagation::config prop_cfg{};
  prop_cfg.navigation.intersection.overstep_tolerance =
      -100.f * unit<float>::um;
  propagator_t p{prop_cfg};

  constexpr scalar q{-1.f};
  constexpr scalar iniP{10.f * unit<scalar>::GeV};

  // Initial track parameters directing x-axis
  bound_parameters_vector<test_algebra> bound_vector{};
  bound_vector.set_theta(constant<scalar>::pi_2);
  bound_vector.set_qop(q / iniP);

  auto bound_cov = matrix::identity<covariance_t>();
  getter::element(bound_cov, e_bound_phi, e_bound_phi) = 0.f;
  getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;

  // bound track parameter
  const bound_track_parameters<test_algebra> bound_param(
      det.surface(0u).identifier(), bound_vector, bound_cov);

  std::size_t n_samples{100000u};
  std::vector<scalar> phis;
  std::vector<scalar> thetas;

  for (std::size_t i = 0u; i < n_samples; i++) {
    pathlimit_aborter_t::state aborter_state{};
    // Seed = sample id
    actor::parameter_updater_state<test_algebra> updater_state{prop_cfg,
                                                               bound_param};
    simulator_t::state simulator_state{i};
    simulator_state.do_energy_loss = false;

    // Create actor states tuples
    auto actor_states =
        detray::tie(aborter_state, updater_state, simulator_state);

    propagator_t::state state(bound_param, det);
    state.debug(true);

    // Propagate the entire detector
    ASSERT_TRUE(p.propagate(state, actor_states));

    const auto& final_param = updater_state.bound_params();

    // Updated phi and theta variance
    if (i == 0u) {
      interactor_t{}.update_angle_variance(
          bound_cov, traj.dir(), simulator_state.projected_scattering_angle);
    }

    phis.push_back(final_param.phi());
    thetas.push_back(final_param.theta());
  }
  scalar phi_variance{statistics::rms(phis, bound_param.phi())};
  scalar theta_variance{statistics::rms(thetas, bound_param.theta())};

  // Get the phi and theta variance
  scalar ref_phi_variance =
      getter::element(bound_cov, e_bound_phi, e_bound_phi);
  scalar ref_theta_variance =
      getter::element(bound_cov, e_bound_theta, e_bound_theta);

  // Tolerate upto 1% difference
  EXPECT_NEAR((phi_variance - ref_phi_variance) / ref_phi_variance, 0.f, 1e-2f);
  EXPECT_NEAR((theta_variance - ref_theta_variance) / ref_theta_variance, 0.f,
              1e-2f);

  // To make sure that the variances are not zero
  EXPECT_TRUE(ref_phi_variance > 1e-9f && ref_theta_variance > 1e-9f);
}

// Material interaction test with telescope Geometry with volume material
GTEST_TEST(detray_material, telescope_geometry_volume_material) {
  vecmem::host_memory_resource host_mr;

  // Propagator types
  using bfield_t = bfield::const_field_t<scalar>;
  using stepper_t = rk_stepper<bfield_t::view_t, test_algebra>;

  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar>;
  using actor_chain_t = actor_chain<pathlimit_aborter_t>;

  constexpr std::size_t cache_size{navigation::default_cache_size};

  // Bfield setup
  test::vector3 B_z{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                    2.f * unit<scalar>::T};
  const bfield_t const_bfield = create_const_field<scalar>(B_z);

  // Track setup
  constexpr scalar q{-1.f};
  constexpr scalar iniP{10.f * unit<scalar>::GeV};

  bound_parameters_vector<test_algebra> bound_vector{};
  bound_vector.set_theta(constant<scalar>::pi_2);
  bound_vector.set_qop(q / iniP);

  auto bound_cov = matrix::identity<covariance_t>();
  getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 0.f;

  // Create actor states tuples
  const scalar path_limit = 100 * unit<scalar>::mm;

  // Build in x-direction from given module positions
  detail::ray<test_algebra> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
  std::vector<scalar> positions = {0.f, 5000.f * unit<scalar>::mm};

  // NO material at modules
  const auto module_mat = vacuum<scalar>();

  // Create telescope geometry
  tel_det_config<test_algebra, rectangle2D> tel_cfg{50000.f * unit<scalar>::mm,
                                                    50000.f * unit<scalar>::mm};
  tel_cfg.positions(positions).pilot_track(traj).module_material(module_mat);

  std::vector<material<scalar>> vol_mats = {
      vacuum<scalar>(), isobutane<scalar>(), silicon<scalar>(),
      tungsten<scalar>()};

  for (const auto& mat : vol_mats) {
    tel_cfg.volume_material(mat);
    const auto [det, names] =
        build_telescope_detector<test_algebra>(host_mr, tel_cfg);

    // bound track parameter at first physical plane
    const bound_track_parameters<test_algebra> bound_param(
        det.surface(0u).identifier(), bound_vector, bound_cov);

    using navigator_t = caching_navigator<decltype(det), cache_size,
                                          navigation::print_inspector>;

    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator is built from the stepper and navigator
    propagation::config prop_cfg{};
    prop_cfg.navigation.intersection.overstep_tolerance =
        -100.f * unit<float>::um;
    propagator_t p{prop_cfg};

    propagator_t::stepper_type::magnetic_field_type bfield_view(const_bfield);
    propagator_t::state state(bound_param, bfield_view, det);

    pathlimit_aborter_t::state abrt_state{path_limit};
    auto actor_states = detray::tie(abrt_state);

    p.propagate(state, actor_states);

    const auto newP = state.stepping()().p(ptc.charge());
    const auto mass = ptc.mass();

    const auto eloss_approx =
        interaction<scalar>().compute_energy_loss_bethe_bloch(
            state.stepping().path_length(), mat, ptc, {ptc, bound_param.qop()});

    const auto iniE = std::sqrt(iniP * iniP + mass * mass);
    const auto newE = std::sqrt(newP * newP + mass * mass);
    const auto eloss = iniE - newE;

    if (mat == vacuum<scalar>()) {
      ASSERT_FLOAT_EQ(float(eloss), 0.f);
    } else {
      ASSERT_TRUE(eloss > 0.f) << "mat.: " << mat << "\n"
                               << state.navigation().inspector().to_string();
    }

    ASSERT_NEAR(eloss, eloss_approx, eloss * 0.01);
  }
}
