/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/actors/ckf_aborter.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"
#include "traccc/fitting/kalman_filter/is_line_visitor.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/sanity/contiguous_on.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <detray/propagator/actors/parameter_updater.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace traccc::host::details {

/// Templated implementation of the track finding algorithm.
///
/// Concrete track finding algorithms can use this function with the appropriate
/// specializations, to find tracks on top of a specific detector type, magnetic
/// field type, and track finding configuration.
///
/// @tparam detector_t The (host) detector type to use
/// @tparam bfield_t   The magnetic field type to use
///
/// @param det               The detector object
/// @param field             The magnetic field object
/// @param measurements_view All measurements in an event
/// @param seeds_view        All seeds in an event to start the track finding
///                          with
/// @param config            The track finding configuration
/// @param mr                The memory resource to use
/// @param log               The logger object to use
///
/// @return A container of the found tracks
///
template <typename detector_t, typename bfield_t>
edm::track_container<typename detector_t::algebra_type>::host
combinatorial_kalman_filter(
    const detector_t& det, const bfield_t& field,
    const edm::measurement_collection::const_view& measurements_view,
    const bound_track_parameters_collection_types::const_view& seeds_view,
    const finding_config& config, vecmem::memory_resource& mr,
    const Logger& /*log*/) {
  TRACCC_VERBOSE_HOST_DEVICE("Running CKF...");

  assert(
      config.min_step_length_for_next_surface >
          math::fabs(
              config.propagation.navigation.intersection.overstep_tolerance) &&
      "Min step length for the next surface should be higher than the "
      "overstep tolerance");
  assert(config.min_track_candidates_per_track >= 1);

  /// The algebra type
  using algebra_type = typename detector_t::algebra_type;
  /// The scalar type
  using scalar_type = detray::dscalar<algebra_type>;

  // Create a logger.
  // @TODO: Turn back on, once detray can use the ACTS logger
  // auto logger = [&log]() -> const Logger& { return log; };

  /*****************************************************************
   * Measurement Operations
   *****************************************************************/

  // Create the measurement container.
  edm::measurement_collection::const_device measurements{measurements_view};

  // Check contiguity of the measurements
  assert(is_contiguous_on([](const auto& value) { return value; },
                          measurements.surface_link()));

  // Get index ranges in the measurement container per detector surface
  std::vector<unsigned int> meas_ranges;
  meas_ranges.reserve(det.surfaces().size());

  for (const typename detector_t::surface_type& sf_desc : det.surfaces()) {
    // Measurements can only be found on sensitive surfaces
    if (!sf_desc.is_sensitive()) {
      // Lower range index is the upper index of the previous range
      // This is guaranteed by the measurement sorting step
      const unsigned int sf_idx{sf_desc.index()};
      const unsigned int lo{sf_idx == 0u ? 0u : meas_ranges[sf_idx - 1u]};

      // Hand the upper index of the previous range through to assign
      // the lower index of the next valid range correctly
      meas_ranges.push_back(lo);
      continue;
    }

    auto up = std::upper_bound(measurements.surface_link().begin(),
                               measurements.surface_link().end(),
                               sf_desc.identifier());
    meas_ranges.push_back(static_cast<unsigned int>(
        std::distance(measurements.surface_link().begin(), up)));
  }

  const edm::measurement_collection::const_device::size_type n_meas =
      measurements.size();

  std::vector<std::vector<candidate_link>> links;
  links.resize(config.max_track_candidates_per_track);

  std::vector<std::vector<std::size_t>> param_to_link;
  param_to_link.resize(config.max_track_candidates_per_track);

  std::vector<std::pair<unsigned int, unsigned int>> tips;

  // Create propagator
  detray::propagation::config prop_cfg{config.propagation};
  prop_cfg.navigation.estimate_scattering_noise = false;
  traccc::details::ckf_propagator_t<detector_t, bfield_t> propagator(prop_cfg);

  // Create the input seeds container.
  bound_track_parameters_collection_types::const_device seeds{seeds_view};

  // Copy seed to input parameters
  std::vector<bound_track_parameters<algebra_type>> in_params(seeds.size());
  std::copy(seeds.begin(), seeds.end(), in_params.begin());
  std::vector<unsigned int> n_trks_per_seed(seeds.size());

  std::vector<bound_track_parameters<algebra_type>> out_params;

  for (unsigned int step = 0u; step < config.max_track_candidates_per_track;
       step++) {
    TRACCC_VERBOSE_HOST("Starting step "
                        << step + 1 << " / "
                        << config.max_track_candidates_per_track);

    // Iterate over input parameters
    const std::size_t n_in_params = in_params.size();
    TRACCC_VERBOSE_HOST("No. Params: " << n_in_params);

    // Terminate if there is no parameter to proceed
    if (n_in_params == 0) {
      break;
    }

    // Rough estimation on out parameters size
    out_params.reserve(n_in_params);

    // Previous step ID
    std::fill(n_trks_per_seed.begin(), n_trks_per_seed.end(), 0u);

    // Parameters updated by Kalman fitter
    std::vector<bound_track_parameters<algebra_type>> updated_params;

    for (unsigned int in_param_id = 0; in_param_id < n_in_params;
         in_param_id++) {
      bound_track_parameters<algebra_type>& in_param = in_params[in_param_id];

      assert(!in_param.is_invalid());

      const unsigned int orig_param_id =
          (step == 0 ? in_param_id
                     : links[step - 1][param_to_link[step - 1][in_param_id]]
                           .seed_idx);
      const unsigned int skip_counter =
          (step == 0 ? 0
                     : links[step - 1][param_to_link[step - 1][in_param_id]]
                           .n_skipped);
      const unsigned int consecutive_skipped =
          (step == 0 ? 0
                     : links[step - 1][param_to_link[step - 1][in_param_id]]
                           .n_consecutive_skipped);
      const scalar prev_chi2_sum =
          (step == 0 ? 0.f
                     : links[step - 1][param_to_link[step - 1][in_param_id]]
                           .chi2_sum);
      const unsigned int prev_ndf_sum =
          (step == 0
               ? 0
               : links[step - 1][param_to_link[step - 1][in_param_id]].ndf_sum);

      TRACCC_VERBOSE_HOST("Processing input parameter "
                          << in_param_id + 1 << " / " << n_in_params);
      TRACCC_VERBOSE_HOST("-> orig_param_id="
                          << orig_param_id << ", skip_counter=" << skip_counter
                          << "\nVec:\n"
                          << in_param.vector());
      TRACCC_DEBUG_HOST("Cov:\n" << in_param.covariance());

      // Get surface corresponding to bound params
      const detray::tracking_surface sf{det, in_param.surface_link()};

      TRACCC_VERBOSE_HOST(" Free params:\n"
                          << sf.bound_to_free_vector({}, in_param));

      /*****************************************************************
       * Find tracks (CKF)
       *****************************************************************/

      // Iterate over the measurements for this surface
      const unsigned int sf_idx{sf.index()};
      const unsigned int lo{sf_idx == 0u ? 0u : meas_ranges[sf_idx - 1]};
      const unsigned int up{meas_ranges[sf_idx]};

      std::vector<
          std::tuple<candidate_link, bound_track_parameters<algebra_type>>>
          best_links;

      const bool is_line = detail::is_line(sf);

      // Iterate over the measurements
      TRACCC_VERBOSE_HOST("No. measurements: " << (up - lo));
      for (unsigned int meas_id = lo; meas_id < up; meas_id++) {
        TRACCC_VERBOSE_HOST("Testing measurement: " << meas_id);

        // The measurement on surface to handle.
        const edm::measurement meas = measurements.at(meas_id);

        const scalar_type chi2 = measurement_selector::predicted_chi2(
            meas, in_param, config.meas_calibration, is_line);

        // If the measurement is outside the chi2 cut, skip it
        if (chi2 > config.chi2_max || chi2 < 0.f) {
          continue;
        }

        // Create a standalone track state object.
        edm::track_state trk_state =
            edm::make_track_state<algebra_type>(measurements, meas_id);
        trk_state.filtered_chi2() = chi2;

        // Kalman filter status code
        kalman_fitter_status res{kalman_fitter_status::ERROR_OTHER};

        // Don't run the filter on the first measurement
        if (step == 0 && !sf.has_material()) {
          // Only do this for the actual seed measurement
          if (chi2 == 0.f) {
            res = kalman_fitter_status::SUCCESS;

            // Copy the full track parameters
            // TODO: Apply calibration ?
            trk_state.filtered_params() = in_param;

            // Update measurement covariance
            const auto V =
                measurement_selector::calibrated_measurement_covariance<
                    algebra_type, 2>(meas, config.meas_calibration);

            auto& filtered_cov = trk_state.filtered_params().covariance();
            getter::element(filtered_cov, e_bound_loc0, e_bound_loc0) =
                getter::element(V, 0, 0);
            getter::element(filtered_cov, e_bound_loc1, e_bound_loc1) =
                getter::element(V, 1, 1);
          }
        } else {
          // Run the Kalman update on the track state
          constexpr gain_matrix_updater<algebra_type> kalman_updater{};
          res = kalman_updater(trk_state, meas, in_param,
                               config.meas_calibration, is_line);
        }

        TRACCC_DEBUG_HOST("KF status: " << fitter_debug_msg{res}());

        // The chi2 from Kalman update should be less than chi2_max
        if (res == kalman_fitter_status::SUCCESS) {
          TRACCC_VERBOSE_HOST("Found measurement: " << meas_id);

          best_links.push_back({{.step = step,
                                 .previous_candidate_idx = in_param_id,
                                 .meas_idx = meas_id,
                                 .seed_idx = orig_param_id,
                                 .n_skipped = skip_counter,
                                 .n_consecutive_skipped = 0,
                                 .chi2 = chi2,
                                 .chi2_sum = prev_chi2_sum + chi2,
                                 .ndf_sum = prev_ndf_sum + meas.dimensions()},
                                trk_state.filtered_params()});
        }
      }

      // Sort the links by chi2
      std::sort(best_links.begin(), best_links.end(),
                [](const auto& a, const auto& b) {
                  return std::get<0>(a).chi2 < std::get<0>(b).chi2;
                });
      // Take the best links
      const unsigned int n_branches =
          std::min(config.max_num_branches_per_surface,
                   static_cast<unsigned int>(best_links.size()));
      TRACCC_VERBOSE_HOST("Found " << n_branches << " branches for step "
                                   << step << " and input parameter "
                                   << in_param_id + 1);
      for (unsigned int i = 0; i < n_branches; ++i) {
        const auto& [link, filtered_params] = best_links[i];

        // Add the link to the links container
        links[step].push_back(link);

        // Add the updated parameter to the updated parameters
        updated_params.push_back(filtered_params);
        TRACCC_DEBUG_HOST("updated_params[" << updated_params.size() - 1
                                            << "] = " << updated_params.back());
      }

      /*****************************************************************
       * Add a dummy links in case of no branches
       *****************************************************************/

      if (n_branches == 0) {
        // Put an invalid link with max item id
        links[step].push_back(
            {.step = step,
             .previous_candidate_idx = in_param_id,
             .meas_idx = std::numeric_limits<unsigned int>::max(),
             .seed_idx = orig_param_id,
             .n_skipped = skip_counter + 1,
             .n_consecutive_skipped = consecutive_skipped + 1,
             .chi2 = std::numeric_limits<traccc::scalar>::max(),
             .chi2_sum = prev_chi2_sum,
             .ndf_sum = prev_ndf_sum});

        TRACCC_VERBOSE_HOST("Hole state created");

        updated_params.push_back(in_param);
        TRACCC_DEBUG_HOST("updated_params[" << updated_params.size() - 1
                                            << "] = " << updated_params.back());
      }
    }

    /*
     * Track deduplication.
     *
     * For documentation, see the device version.
     */
    const std::size_t n_links = links[step].size();
    std::vector<unsigned int> param_liveness;
    param_liveness.resize(n_links);

    for (std::size_t i = 0; i < param_liveness.size(); ++i) {
      param_liveness.at(i) = 1u;
    }

    if (step >= config.duplicate_removal_minimum_length) {
      std::map<std::size_t, std::vector<std::size_t>> last_meas_to_tracks_map;

      for (std::size_t i = 0; i < n_links; ++i) {
        candidate_link L = links.at(step).at(i);

        while (L.meas_idx >= n_meas && L.step != 0u) {
          const std::size_t link_pos =
              param_to_link.at(L.step - 1u).at(L.previous_candidate_idx);
          L = links.at(L.step - 1u).at(link_pos);
        }

        last_meas_to_tracks_map[L.meas_idx].push_back(i);
      }

      for (const auto& [_, tracks] : last_meas_to_tracks_map) {
        for (std::size_t i = 0; i < tracks.size(); ++i) {
          const candidate_link& Lthisbase = links.at(step).at(tracks.at(i));
          const scalar prob_this = prob(
              Lthisbase.chi2_sum, static_cast<scalar>(Lthisbase.ndf_sum - 5));

          if (step + 1 - Lthisbase.n_skipped <=
                  config.duplicate_removal_minimum_length ||
              Lthisbase.ndf_sum <= 5) {
            continue;
          }

          for (std::size_t j = 0; j < tracks.size(); ++j) {
            if (i == j) {
              continue;
            }

            candidate_link Lthis = Lthisbase;
            candidate_link Lthat = links.at(step).at(tracks.at(j));

            if (step + 1 - Lthat.n_skipped <=
                    config.duplicate_removal_minimum_length ||
                Lthisbase.ndf_sum <= 5) {
              continue;
            }

            bool this_is_dominated = true;

            const scalar prob_that =
                prob(Lthat.chi2_sum, static_cast<scalar>(Lthat.ndf_sum - 5));

            while (true) {
              while (Lthis.meas_idx >= n_meas && Lthis.step != 0u) {
                const std::size_t link_pos =
                    param_to_link.at(Lthis.step - 1u)
                        .at(Lthis.previous_candidate_idx);

                Lthis = links.at(Lthis.step - 1u).at(link_pos);
              }
              while (Lthat.meas_idx >= n_meas && Lthat.step != 0u) {
                const std::size_t link_pos =
                    param_to_link.at(Lthat.step - 1u)
                        .at(Lthat.previous_candidate_idx);

                Lthat = links.at(Lthat.step - 1u).at(link_pos);
              }

              if (Lthis.meas_idx == Lthat.meas_idx) {
                if (Lthis.step == 0) {
                  break;
                } else if (Lthat.step == 0) {
                  this_is_dominated = false;
                  break;
                } else {
                  const std::size_t link_pos_this =
                      param_to_link.at(Lthis.step - 1u)
                          .at(Lthis.previous_candidate_idx);
                  Lthis = links.at(Lthis.step - 1u).at(link_pos_this);
                  const std::size_t link_pos_that =
                      param_to_link.at(Lthat.step - 1u)
                          .at(Lthat.previous_candidate_idx);
                  Lthat = links.at(Lthat.step - 1u).at(link_pos_that);
                }
              } else {
                this_is_dominated = false;
                break;
              }
            }

            if (prob_this != prob_that) {
              this_is_dominated &= prob_that >= prob_this;
            } else {
              this_is_dominated &= tracks.at(j) < tracks.at(i);
            }

            if (this_is_dominated) {
              TRACCC_VERBOSE_HOST("Track is dead (deduplication)!");
              param_liveness.at(tracks.at(i)) = 0u;
              break;
            }
          }
        }
      }
    }

    /*********************************
     * Propagate to the next surface
     *********************************/
    for (unsigned int link_id = 0; link_id < n_links; link_id++) {
      if (param_liveness.at(link_id) == 0u) {
        continue;
      }

      const unsigned int seed_idx = links.at(step).at(link_id).seed_idx;
      n_trks_per_seed[seed_idx]++;

      if (n_trks_per_seed[seed_idx] > config.max_num_branches_per_seed) {
        continue;
      }

      // If number of skips is larger than the maximum value, consider the
      // link to be a tip
      if (links.at(step).at(link_id).n_skipped >
          config.max_num_skipping_per_cand) {
        TRACCC_WARNING_HOST(
            "Create tip: Max no. of holes reached! Bound param:\n"
            << updated_params[link_id].vector());
        tips.push_back({step, link_id});
        continue;
      }

      // If number of consecutive skips is larger than the maximum value,
      // consider the link to be a tip
      if (links.at(step).at(link_id).n_consecutive_skipped >
          config.max_num_consecutive_skipped) {
        TRACCC_WARNING_HOST(
            "Create tip: Max no. of consecutive holes reached! Bound "
            "param:\n"
            << updated_params[link_id].vector());
        tips.push_back({step, link_id});
        continue;
      }

      const bound_track_parameters<algebra_type>& param =
          updated_params[link_id];

      // Create propagator state
      typename traccc::details::ckf_propagator_t<detector_t, bfield_t>::state
          propagation(param, field, det);
      propagation.set_particle(
          detail::correct_particle_hypothesis(config.ptc_hypothesis, param));

      propagation.stepping()
          .template set_constraint<detray::step::constraint::e_accuracy>(
              config.propagation.stepping.step_constraint);

      typename detray::actor::pathlimit_aborter<scalar_type>::state
          aborter_state;
      detray::actor::parameter_updater_state<typename detector_t::algebra_type>
          updater_state{prop_cfg, param};
      traccc::details::ckf_interactor_t::state interactor_state;
      typename detray::actor::momentum_aborter<scalar_type>::state
          momentum_aborter_state{};
      typename ckf_aborter::state ckf_aborter_state;

      // Update the actor config
      // Notify the KF and material interaction only at the first
      // propagation initialization. For all subsequent initializations,
      // they will have run already on the previous step
      updater_state.notify_on_initial(step == 0);
      momentum_aborter_state.min_pT(static_cast<scalar_type>(config.min_pT));
      momentum_aborter_state.min_p(static_cast<scalar_type>(config.min_p));
      ckf_aborter_state.min_step_length =
          config.min_step_length_for_next_surface;
      ckf_aborter_state.max_count = config.max_step_counts_for_next_surface;

      // Propagate to the next surface
      TRACCC_DEBUG_HOST("Propagating... ");
      propagator.propagate(
          propagation,
          detray::tie(aborter_state, updater_state, interactor_state,
                      momentum_aborter_state, ckf_aborter_state));
      TRACCC_DEBUG_HOST("Finished propagation");

      // If a surface found, add the parameter for the next step
      bool valid_track{ckf_aborter_state.success};
      if (valid_track) {
        assert(propagation.navigation().is_on_sensitive());
        assert(!updater_state.bound_params().is_invalid());
        TRACCC_DEBUG_HOST(
            "On surface: " << propagation.navigation().geometry_identifier());

        const bound_track_parameters<algebra_type>& out_param =
            updater_state.bound_params();

        const scalar theta = out_param.theta();
        if (theta <= 0.f || theta >= 2.f * constant<traccc::scalar>::pi) {
          TRACCC_ERROR_HOST("Theta is hit pole after propagation");
          valid_track = false;
        }

        if (!std::isfinite(out_param.phi())) {
          TRACCC_ERROR_HOST(
              "Phi is infinite after propagation (Matrix inversion)");
          valid_track = false;
        }

        if (math::fabs(out_param.qop()) == 0.f) {
          TRACCC_ERROR_HOST("q over p is zero after propagation");
          valid_track = false;
        }

        if (valid_track) {
          out_params.push_back(out_param);
          param_to_link[step].push_back(link_id);
        }
      }
      // Unless the track found a surface, it is considered a
      // tip
      if (!valid_track &&
          (step >= (config.min_track_candidates_per_track - 1u))) {
        if (!ckf_aborter_state.success) {
          // HACK: Silence SonarCloud S3923
          static_cast<void>(ckf_aborter_state.success);
          TRACCC_VERBOSE_HOST("Create tip: No next sensitive found");
        } else {
          TRACCC_VERBOSE_HOST("Create tip: Encountered error");
        }
        tips.push_back({step, link_id});
      }

      // If no more CKF step is expected, current candidate is
      // kept as a tip
      if (ckf_aborter_state.success &&
          (step == (config.max_track_candidates_per_track - 1u))) {
        TRACCC_ERROR_HOST("Create tip: Max no. candidates");
        tips.push_back({step, link_id});
      }
    }

    in_params = std::move(out_params);
    out_params.clear();
  }

  /**********************
   * Build tracks
   **********************/

  // Number of found tracks = number of tips
  typename edm::track_container<algebra_type>::host output_candidates{
      mr, measurements_view};
  output_candidates.tracks.reserve(tips.size());

  for (const auto& tip : tips) {
    // Get the link corresponding to tip
    candidate_link L = links.at(tip.first).at(tip.second);

    const unsigned int n_cands = tip.first + 1 - L.n_skipped;

    // Skip if the number of tracks candidates is too small
    if (n_cands < config.min_track_candidates_per_track ||
        n_cands > config.max_track_candidates_per_track) {
      continue;
    }

    vecmem::vector<edm::track_constituent_link> cands_per_track;
    cands_per_track.resize(n_cands);

    // Track summary variables
    scalar ndf_sum = 0.f;
    scalar chi2_sum = 0.f;

    // Reversely iterate to fill the track candidates
    for (auto it = cands_per_track.rbegin(); it != cands_per_track.rend();
         it++) {
      while (L.meas_idx >= n_meas && L.step != 0u) {
        const std::size_t link_pos =
            param_to_link.at(L.step - 1u).at(L.previous_candidate_idx);

        L = links.at(L.step - 1u).at(link_pos);
      }

      // Break if the measurement is still invalid
      if (L.meas_idx >= measurements.size()) {
        break;
      }

      *it = {edm::track_constituent_link::measurement, L.meas_idx};

      // Sanity check on chi2
      assert(L.chi2 < std::numeric_limits<traccc::scalar>::max());
      assert(L.chi2 >= 0.f);

      ndf_sum += static_cast<scalar>(measurements.at(it->index).dimensions());
      chi2_sum += L.chi2;

      // Break the loop if the iterator is at the first candidate and
      // fill the seed
      if (it == cands_per_track.rend() - 1) {
        const bound_track_parameters<algebra_type>& cand_seed =
            seeds.at(L.seed_idx);
        ndf_sum = ndf_sum - 5.f;
        const scalar_type pval = prob(chi2_sum, ndf_sum);

        // Add seed and track candidates to the output container
        output_candidates.tracks.push_back({});
        edm::track track =
            output_candidates.tracks.at(output_candidates.tracks.size() - 1);
        track.fit_outcome() = track_fit_outcome::UNKNOWN;
        track.params() = cand_seed;
        track.ndf() = ndf_sum;
        track.chi2() = chi2_sum;
        track.pval() = pval;
        track.nholes() = L.n_skipped;
        track.constituent_links() = cands_per_track;

      } else {
        const std::size_t l_pos =
            param_to_link.at(L.step - 1u).at(L.previous_candidate_idx);

        L = links.at(L.step - 1u).at(l_pos);
      }
    }
  }

  return output_candidates;
}

}  // namespace traccc::host::details
