/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/actors/measurement_kalman_updater.hpp"
#include "traccc/finding/details/progressive_kalman_filter_types.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/finding/track_state_candidate.hpp"
#include "traccc/sanity/contiguous_on.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace traccc::details {

/// Templated implementation of a progressive Kalman Filter for track finding.
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
/// @param config            The track finding configuration
///
/// @return A struct that contains information about the found track
///
template <typename detector_t, typename bfield_t>
TRACCC_HOST_DEVICE inline track_stats<typename detector_t::scalar_type>
progressive_kalman_filter(
    const detector_t& det, const bfield_t& field,
    const typename edm::measurement_collection::const_view& measurements_view,
    const vecmem::data::vector_view<unsigned int> measurement_ranges_view,
    const bound_track_parameters<typename detector_t::algebra_type>& seed,
    [[maybe_unused]] const unsigned int seed_idx,
    void* track_state_candidate_ptr,
    vecmem::data::jagged_vector_view<typename detector_t::surface_type>
        surfaces_view,
    const finding_config& cfg) {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = detray::dscalar<algebra_t>;

    // Detray propagation types
    using propagator_t =
        traccc::details::pkf_propagator_t<detector_t, bfield_t>;

    // Create the measurement container.
    typename edm::measurement_collection::const_device measurements{
        measurements_view};

    // Create detray propagator
    auto prop_cfg{cfg.propagation};
    prop_cfg.navigation.estimate_scattering_noise = false;
    propagator_t propagator(prop_cfg);

    TRACCC_INFO_HOST("Seed: " << seed_idx);
    // Add the information also to the clog stream
    TRACCC_VERBOSE_HOST_DEVICE("Seed: %d", seed_idx);

    const auto ptc_hypothesis{
        detail::correct_particle_hypothesis(cfg.ptc_hypothesis, seed)};

    // Check if the seed should be forwarded to the KF
    const scalar_t q{ptc_hypothesis.charge()};
    if (seed.pT(q) <= static_cast<scalar_t>(cfg.min_pT) / 10.f) {
        TRACCC_WARNING_HOST_DEVICE(
            "Seed below min. transverse momentum: |pT| = %f MeV", seed.pT(q));
        return {};
    }
    if (seed.p(q) <= static_cast<scalar_t>(cfg.min_p) / 10.f) {
        TRACCC_WARNING_HOST_DEVICE("Seed below min. momentum: |p| = %f MeV",
                                   seed.p(q));
        return {};
    }

    // Construct propagation state
    typename propagator_t::state propagation(seed, field, det);
    propagation.set_particle(
        detail::correct_particle_hypothesis(cfg.ptc_hypothesis, seed));

    // Pathlimit aborter
    typename detray::actor::pathlimit_aborter<scalar_t>::state
        path_aborter_state;
    // Momentum aborter
    typename detray::actor::momentum_aborter<scalar_t>::state
        momentum_aborter_state{};
    // Track parameter transporter
    typename detray::actor::parameter_updater_state<algebra_t> updater_state{
        prop_cfg, seed};
    // Material interactor
    typename detray::actor::pointwise_material_interactor<algebra_t>::state
        interactor_state;
    // Do the measurement selection
    typename traccc::measurement_updater<algebra_t>::state meas_updater_state{
        measurements, measurement_ranges_view, track_state_candidate_ptr,
        cfg.run_smoother};
    // Collect the surface geometry identifiers for the Kalman smoother
    typename detray::actor::surface_sequencer<
        typename detector_t::surface_type>::state sequencer_state{
        vecmem::device_vector<typename detector_t::surface_type>(
            *(surfaces_view.ptr() +
              (cfg.run_smoother == smoother_type::e_kalman ? seed_idx : 0u)))};

    path_aborter_state.set_path_limit(cfg.propagation.stepping.path_limit);
    momentum_aborter_state.min_pT(static_cast<scalar_t>(cfg.min_pT));
    momentum_aborter_state.min_p(static_cast<scalar_t>(cfg.min_p));

    meas_updater_state.max_chi2 = cfg.chi2_max;
    meas_updater_state.max_n_track_states =
        static_cast<std::uint_least16_t>(cfg.max_track_candidates_per_track);
    meas_updater_state.max_n_holes =
        static_cast<std::uint_least16_t>(cfg.max_num_skipping_per_cand);
    meas_updater_state.max_n_consecutive_holes =
        static_cast<std::uint_least16_t>(cfg.max_num_consecutive_skipped);
    meas_updater_state.m_calib_cfg = cfg.meas_calibration;

    auto actor_states = detray::tie(path_aborter_state, sequencer_state,
                                    updater_state, interactor_state,
                                    meas_updater_state, momentum_aborter_state);

    assert(meas_updater_state.m_stats.n_holes < cfg.max_num_skipping_per_cand);
    assert(meas_updater_state.m_stats.n_consecutive_holes <
           cfg.max_num_consecutive_skipped);

    // Do not rerun measurement selection when resuming
    updater_state.notify_on_initial(true);
    propagator.propagate(propagation, actor_states);

    bool good_track = propagator.finished(propagation);

    // Check if the track should be continued
    const auto& free_param = propagation.stepping()();
    const scalar_t q_new{propagation.stepping().particle_hypothesis().charge()};
    if (good_track &&
        free_param.pT(q_new) <= static_cast<scalar_t>(cfg.min_pT)) {
        TRACCC_WARNING_HOST_DEVICE(
            "Track below min. transverse momentum: |pT| = %f MeV",
            free_param.pT(q_new));
        return {};
    }
    if (good_track && free_param.p(q_new) <= static_cast<scalar_t>(cfg.min_p)) {
        TRACCC_WARNING_HOST_DEVICE("Track below min. momentum: |p| = %f MeV",
                                   free_param.p(q_new));
        return {};
    }

    TRACCC_VERBOSE_HOST_DEVICE("Track following finished. Building track...");

    const traccc::track_stats<scalar_t>& trk_stats = meas_updater_state.m_stats;
    const unsigned int n_track_states{trk_stats.n_track_states};

    TRACCC_DEBUG_HOST("Seed params:\n" << seed);
    TRACCC_VERBOSE_HOST_DEVICE(
        "Found track for seed %d: (%d track states, %d holes)", seed_idx,
        n_track_states, trk_stats.n_holes);

    assert(n_track_states <= cfg.max_track_candidates_per_track);
    if (good_track && n_track_states < cfg.min_track_candidates_per_track) {
        TRACCC_WARNING_HOST_DEVICE("Short track (%d track states): discarding",
                                   n_track_states);
        good_track = false;
    }

    // Check track stats
    const int ndf_sum{static_cast<int>(trk_stats.ndf_sum) - 5};
    if (good_track && ndf_sum < 0) {
        TRACCC_ERROR_HOST_DEVICE("Negative NDF sum for track: discarding");
        good_track = false;
    }

    if (good_track && trk_stats.n_holes > cfg.max_num_skipping_per_cand + 1) {
        TRACCC_ERROR_HOST_DEVICE("Hole count incorrect: max = %d, found = %d",
                                 cfg.max_num_skipping_per_cand,
                                 trk_stats.n_holes);
    }

    if (good_track &&
        trk_stats.n_consecutive_holes > cfg.max_num_consecutive_skipped + 1) {
        TRACCC_ERROR_HOST_DEVICE(
            "Consecutive hole count incorrect: max = %d, found = %d",
            cfg.max_num_consecutive_skipped, trk_stats.n_consecutive_holes);
    }

    // Return the collected track statistics for further evaluation
    return good_track ? trk_stats : traccc::track_stats<scalar_t>{};
}

}  // namespace traccc::details
