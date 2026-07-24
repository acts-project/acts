// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/utils/perigee_stopper.hpp"
#include "detray/test/validation/material_validation_utils.hpp"
#include "detray/test/validation/navigation_validation_utils.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <iostream>
#include <memory>
#include <optional>
#include <string>

namespace detray {

template <typename detector_t>
using intersection_type =
    intersection2D<typename detector_t::surface_type,
                   typename detector_t::algebra_type, true>;

template <typename detector_t>
using candidate_type =
    navigation::detail::candidate_record<intersection_type<detector_t>>;

/// Configure the Kalman Filter comparison test
struct propagation_validation_config {
  propagation::config propagation;

  pdg_particle<float> particle = detray::muon<float>();

  bool do_multiple_scattering{true};
  bool do_energy_loss{true};

  bool display_svg{true};
  float max_percent_missed{1.f};
  float max_percent_additional{15.f};
};

/// Run a detray propagation and compare the surfaces that were encountered
/// against a reference of truth tracks givben as intersection traces
///
/// @param det the detector
/// @param names the detector and volume names
/// @param bfield the magnetic field representation
/// @param cfg the full test configuration
/// @param tracks the initial particle vertices as free track parameters
/// @param truth_traces_fw the reference data for each encountered module
///
/// @returns whether the validation was successful
template <typename detector_t, typename bfield_view_t>
bool propagation_validation(
    const detector_t& det, const typename detector_t::name_map& names,
    const std::optional<bfield_view_t>& bfield,
    const propagation_validation_config& cfg,
    const std::vector<free_track_parameters<typename detector_t::algebra_type>>&
        tracks,
    std::vector<vecmem::vector<candidate_type<detector_t>>>& truth_traces_fw) {
  // 'false' if any failures were detected
  bool test_successful{true};

  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using stepper_t =
      rk_stepper<bfield_view_t, algebra_t, constrained_step<scalar_t>,
                 stepper_rk_policy<scalar_t>, stepping::print_inspector>;
  using sf_candidate_t =
      traccc::propagation_validator::candidate_type<detector_t>;

  // Host memory resource
  vecmem::host_memory_resource host_mr;

  // Geometry context
  const typename detector_t::geometry_context ctx{};

  // Configure the test
  test::navigation_validation_config<algebra_t> test_cfg{};
  test_cfg.n_tracks(tracks.size());
  test_cfg.ptc_hypothesis(cfg.particle);
  test_cfg.collect_sensitives_only(true);
  test_cfg.fail_on_diff(false);
  test_cfg.display_svg(cfg.display_svg);
  test_cfg.display_only_missed(true);
  test_cfg.verbose(false);

  // Create the backward truth traces by reverting the order
  std::vector<vecmem::vector<sf_candidate_t>> truth_traces_bw;
  truth_traces_bw.reserve(tracks.capacity());

  for (const auto& truth_trace_fw : truth_traces_fw) {
    // Revert the forward trace for the backward propagation
    vecmem::vector<sf_candidate_t> truth_trace_bw(truth_trace_fw.size());
    std::ranges::reverse_copy(truth_trace_fw, truth_trace_bw.begin());

    assert(!truth_trace_bw.empty());

    truth_traces_bw.push_back(std::move(truth_trace_bw));
  }

  // Check truth data
  if (truth_traces_fw.empty() || truth_traces_bw.empty()) {
    DETRAY_ERROR_HOST("Propagation truth data empty");
    return false;
  }

  // Make a tuple of references from a tuple
  auto setup_actor_states = []<typename... T>(detray::dtuple<T...>& t) {
    return detray::tie(detray::detail::get<T>(t)...);
  };

  // Define the actors
  using interactor = actor::pointwise_material_interactor<algebra_t>;
  using perigee_stopper = perigee_stopper<algebra_t>;

  // Reusable actor states
  typename interactor::state interactor_state{};
  interactor_state.do_multiple_scattering = cfg.do_multiple_scattering;
  interactor_state.do_energy_loss = cfg.do_energy_loss;

  typename perigee_stopper::state stopper_state{};

  std::cout << "-----------------------------------"
            << "\nFORWARD - No KF" << std::endl
            << "-----------------------------------\n";

  using parameter_updater = actor::parameter_updater<algebra_t, interactor>;

  using actor_chain_t = actor_chain<parameter_updater>;

  typename parameter_updater::state updater_state{};
  updater_state.noise_estimation_cfg().n_stddev =
      cfg.propagation.navigation.n_scattering_stddev;
  updater_state.noise_estimation_cfg().accumulated_error =
      cfg.propagation.navigation.accumulated_error;
  updater_state.noise_estimation_cfg().estimate_scattering_noise =
      cfg.propagation.navigation.estimate_scattering_noise;

  // Prepare the actor states for every track
  vecmem::vector<typename actor_chain_t::state_tuple> state_tuples{};
  vecmem::vector<typename actor_chain_t::state_ref_tuple> state_ref_tuples{};
  state_tuples.reserve(tracks.size());
  state_ref_tuples.reserve(tracks.size());

  for (std::size_t i = 0u; i < tracks.size(); ++i) {
    // Define the initial covariance
    updater_state.init(tracks.at(i));

    // Copy the configured updater state for this track
    state_tuples.push_back(detray::make_tuple(updater_state, interactor_state));
    state_ref_tuples.push_back(setup_actor_states(state_tuples.back()));
  }

  // Forward navigation
  test_cfg.name(det.name(names) + "_GeV_fw");
  test_cfg.navigation_direction(navigation::direction::e_forward);
  const auto [trk_stats_fw, n_surfaces_fw, n_miss_nav_fw, n_miss_truth_fw,
              step_traces_fw, mat_traces_fw, mat_records_fw] =
      navigation_validator::compare_to_navigation<stepper_t, parameter_updater>(
          test_cfg, host_mr, det, names, ctx, bfield, cfg.propagation,
          truth_traces_fw, tracks, state_ref_tuples);

  // Check stats
  auto n_tracks{static_cast<double>(trk_stats_fw.n_tracks)};
  if (static_cast<double>(trk_stats_fw.n_tracks_w_holes) / n_tracks >
      cfg.max_percent_missed / 100.) {
    DETRAY_ERROR_HOST("Too many tracks with missing surfaces");
    test_successful = false;
  }
  if (static_cast<double>(trk_stats_fw.n_tracks_w_extra) / n_tracks >
      cfg.max_percent_additional / 100.) {
    DETRAY_ERROR_HOST("Too many tracks with additional surfaces");
    test_successful = false;
  }

  std::cout << "\n-----------------------------------\n"
            << "BACKWARD - No KF" << std::endl
            << "-----------------------------------\n";

  using actor_chain_bw_t = actor_chain<parameter_updater, perigee_stopper>;

  // Prepare the actor states for every track
  vecmem::vector<typename actor_chain_bw_t::state_tuple> bw_state_tuples{};
  vecmem::vector<typename actor_chain_bw_t::state_ref_tuple>
      bw_state_ref_tuples{};
  bw_state_tuples.reserve(tracks.size());
  bw_state_ref_tuples.reserve(tracks.size());

  for (std::size_t i = 0u; i < tracks.size(); ++i) {
    // Define the initial covariance
    updater_state.init(tracks.at(i));

    // Copy the configured updater state for this track
    bw_state_tuples.push_back(
        detray::make_tuple(updater_state, interactor_state, stopper_state));
    bw_state_ref_tuples.push_back(setup_actor_states(bw_state_tuples.back()));
  }

  // Backward navigation
  test_cfg.name(det.name(names) + "_GeV_bw");
  test_cfg.navigation_direction(navigation::direction::e_backward);
  const auto [trk_stats_bw, n_surfaces_bw, n_miss_nav_bw, n_miss_truth_bw,
              step_traces_bw, mat_traces_bw, mat_records_bw] =
      navigation_validator::compare_to_navigation<stepper_t, parameter_updater,
                                                  perigee_stopper>(
          test_cfg, host_mr, det, names, ctx, bfield, cfg.propagation,
          truth_traces_bw, tracks, bw_state_ref_tuples);

  // Make sure some data was collected
  assert(trk_stats_fw.n_tracks > 0u);
  assert(n_surfaces_fw.n_total() > 0u);
  assert(trk_stats_bw.n_tracks > 0u);
  assert(n_surfaces_bw.n_total() > 0u);

  assert(trk_stats_fw.n_tracks == trk_stats_bw.n_tracks);

  // Check, the amount of collected material between forward and backward
  assert(step_traces_fw.size() == trk_stats_fw.n_tracks);
  assert(mat_traces_fw.size() == trk_stats_fw.n_tracks);
  assert(mat_records_fw.size() == trk_stats_fw.n_tracks);
  assert(mat_records_fw.size() == mat_records_bw.size());
  assert(step_traces_fw.size() == step_traces_bw.size());
  assert(mat_traces_fw.size() == mat_traces_bw.size());

  // Check stats
  n_tracks = static_cast<double>(trk_stats_bw.n_tracks);
  if (static_cast<double>(trk_stats_bw.n_tracks_w_holes) / n_tracks >
      cfg.max_percent_missed / 100.) {
    DETRAY_ERROR_HOST("Too many tracks with missing surfaces");
    test_successful = false;
  }
  if (static_cast<double>(trk_stats_bw.n_tracks_w_extra) / n_tracks >
      cfg.max_percent_additional / 100.) {
    DETRAY_ERROR_HOST("Too many tracks with additional surfaces");
    test_successful = false;
  }

  std::cout << "\n-----------------------------------\n"
            << "MATERIAL TRACE - No KF" << std::endl
            << "-----------------------------------\n";

  constexpr double rel_mat_error{0.01};

  // Material traces contain different surfaces
  std::size_t n_bad_comp{0u};
  // Overall integrated material differs while surface seq. is identical
  std::size_t n_diff_mat{0u};

  // Loop over tracks
  for (std::size_t i = 0u; i < mat_records_fw.size(); ++i) {
    // No material on that track (e.g. detector model without material)
    if (mat_traces_fw[i].empty() && mat_traces_bw[i].empty()) {
      continue;
    }

    std::remove_cvref_t<decltype(mat_traces_bw[i])> inv_mat_trace_bw{};
    if (!mat_traces_bw[i].empty()) {
      inv_mat_trace_bw.resize(mat_traces_bw[i].size());

      // Revert the backward trace to compare to the forward trace
      std::ranges::reverse_copy(mat_traces_bw[i], inv_mat_trace_bw.begin());

      assert(mat_traces_bw[i].size() == inv_mat_trace_bw.size());
    }

    // Compare the material traces and total integrated material per trk
    const auto [is_bad_comp, is_diff_mat] = material_validator::compare_traces(
        mat_traces_fw[i], mat_records_fw[i], inv_mat_trace_bw,
        mat_records_bw[i], i, rel_mat_error, test_cfg.verbose());

    if (is_bad_comp) {
      n_bad_comp++;
    }
    if (is_diff_mat) {
      n_diff_mat++;
    }
  }

  std::cout << "Total no. tracks with diff. material: "
            << (n_bad_comp + n_diff_mat) << " ("
            << 100. * static_cast<double>(n_bad_comp + n_diff_mat) / n_tracks
            << "%)" << std::endl;

  std::cout << "No. identical tracks with diff. material: " << n_diff_mat
            << " (" << 100. * static_cast<double>(n_diff_mat) / n_tracks << "%)"
            << std::endl;
  std::cout << "-----------------------------------\n" << std::endl;

  // Trigger test failures
  if (n_diff_mat != 0) {
    DETRAY_ERROR_HOST("" << n_diff_mat << " tracks have differing material");
    test_successful = false;
  }

  return test_successful;
}

}  // namespace detray
