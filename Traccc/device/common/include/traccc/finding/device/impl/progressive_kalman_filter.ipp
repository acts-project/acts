/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/details/progressive_kalman_filter.hpp"

#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_constituent_link.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/actors/measurement_kalman_updater.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/finding/track_state_candidate.hpp"
#include "traccc/sanity/contiguous_on.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/propagation.hpp"

// Detray include(s).
#include <detray/utils/tuple_helpers.hpp>

namespace traccc::device {

template <typename propagator_t>
TRACCC_HOST_DEVICE inline void progressive_kalman_filter(
    const global_index_t globalIndex, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const typename propagator_t::stepper_type::magnetic_field_type& field_data,
    vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>
        surfaces_view,
    const progressive_kalman_filter_payload& payload) {
  using detector_t = typename propagator_t::detector_type;
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = detray::dscalar<algebra_t>;

  if (globalIndex >= payload.seeds_view.size()) {
    return;
  }

  // Detector
  detector_t det(det_data);

  // Access to initial bound track parameters
  bound_track_parameters_collection_types::device seeds(payload.seeds_view);
  const bound_track_parameters<algebra_t>& seed = seeds.at(globalIndex);

  // Access to measurements and index ranges per surface
  typename edm::measurement_collection::const_device measurements(
      payload.measurements_view);

  // Set the data pointer to the beginning of the range of the track
  const auto cand_offset{static_cast<unsigned int>(
      globalIndex * cfg.max_track_candidates_per_track)};

  track_state_candidate_data<algebra_t> candidate_data(
      cfg.run_smoother, cand_offset, payload.track_cand_view,
      payload.filtered_track_cand_view, payload.full_track_cand_view);

  // Output tracks and track state collection
  typename edm::track_collection<algebra_t>::device track_candidates(
      payload.tracks_view.tracks);
  typename edm::track_state_collection<algebra_t>::device track_states(
      payload.tracks_view.states);

  assert(globalIndex < track_candidates.size());

  // Run the progressive filter for this seed
  const traccc::track_stats<scalar_t> trk_stats =
      traccc::details::progressive_kalman_filter(
          det, field_data, payload.measurements_view,
          payload.measurement_ranges_view, seed, globalIndex,
          candidate_data.ptr(), surfaces_view, cfg);

  // Check track stats and build the new track object
  const unsigned int n_track_states{trk_stats.n_track_states};
  const int ndf_sum{static_cast<int>(trk_stats.ndf_sum) - 5};
  bool good_track{n_track_states > 0u};

  // Link the new states to the final track
  edm::track track = track_candidates.at(globalIndex);

  track.params() = seed;
  track.ndf() = static_cast<scalar_t>(ndf_sum);
  track.chi2() = trk_stats.chi2_sum;
  track.pval() = prob(trk_stats.chi2_sum, static_cast<scalar_t>(ndf_sum));
  track.nholes() = static_cast<unsigned int>(trk_stats.n_holes);
  track.constituent_links().resize(n_track_states);

  if (good_track) {
    track.fit_outcome() = (cfg.run_smoother != smoother_type::e_none)
                              ? track_fit_outcome::SUCCESS
                              : track_fit_outcome::UNKNOWN;

    const auto track_state_offset{globalIndex *
                                  cfg.max_track_candidates_per_track};

    for (unsigned int state_idx = track_state_offset;
         state_idx < track_state_offset + n_track_states; state_idx++) {
      const unsigned int link_idx{state_idx - track_state_offset};

      TRACCC_DEBUG_DEVICE("Adding track state (local idx %d, global idx %d)",
                          link_idx, state_idx);

      // Intermediate type required to build a view
      traccc::track_state_from_candidate<algebra_t>(
          candidate_data.ptr(), cfg.run_smoother, link_idx, measurements, track,
          payload.tracks_view.states);
    }

    TRACCC_INFO_DEVICE(
        "Added track %d to track container (offset %d, #states %d, #tracks "
        "%d)",
        globalIndex, track_state_offset, n_track_states,
        track_candidates.size());
  } else {
    track.fit_outcome() = track_fit_outcome::FAILURE_FITTER;
  }
}

}  // namespace traccc::device
