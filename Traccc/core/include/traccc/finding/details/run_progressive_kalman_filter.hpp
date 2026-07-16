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
#include "traccc/finding/details/progressive_kalman_filter.hpp"
#include "traccc/finding/details/progressive_kalman_filter_types.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/finding/track_state_candidate.hpp"
#include "traccc/fitting/details/kalman_fitting_types.hpp"
#include "traccc/fitting/kalman_smoother.hpp"
#include "traccc/sanity/contiguous_on.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace traccc::host::details {

/// Run the progressive Kalman filter on the host
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
run_progressive_kalman_filter(
    const detector_t& det, const bfield_t& field,
    const typename edm::measurement_collection::const_view& measurements_view,
    const bound_track_parameters_collection_types::const_view& seeds_view,
    const finding_config& cfg, vecmem::memory_resource& mr,
    const Logger& /*log*/) {
  assert(cfg.min_track_candidates_per_track >= 1);

  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = detray::dscalar<algebra_t>;

  // Create the measurement container.
  typename edm::measurement_collection::const_device measurements{
      measurements_view};

  // Check contiguity of the measurements
  assert(is_contiguous_on([](const auto& value) { return value; },
                          measurements.surface_link()));

  TRACCC_INFO_HOST("Run Track Finding: Progressive Kalman Filter");

  // Get index ranges in the measurement container per detector surface
  std::vector<unsigned int> meas_ranges;
  meas_ranges.reserve(det.surfaces().size());

  for (const auto& sf_desc : det.surfaces()) {
    // Measurements can only be found on sensitive surfaces
    if (!sf_desc.is_sensitive()) {
      // Lower range index is the upper index of the previous range
      // This is guaranteed by the measurement sorting step
      const auto sf_idx{sf_desc.index()};
      const unsigned int lo{sf_idx == 0u ? 0u : meas_ranges[sf_idx - 1u]};

      // Hand the upper index of the previous range through to assign
      // the lower index of the next valid range correctly
      meas_ranges.push_back(lo);
      continue;
    }

    auto up = std::ranges::upper_bound(measurements.surface_link(),
                                       sf_desc.identifier());
    meas_ranges.push_back(static_cast<unsigned int>(
        std::distance(measurements.surface_link().begin(), up)));
  }

  // Create the output track container
  typename edm::track_container<algebra_t>::host track_container{
      mr, measurements_view};

  // Create the input seeds container.
  bound_track_parameters_collection_types::const_device seeds{seeds_view};
  const unsigned int n_seeds{seeds.size()};

  if (n_seeds == 0u) {
    TRACCC_DEBUG_HOST("No input seeds: Quitting");
    return track_container;
  }

  // Number of found tracks = number of seeds
  track_container.tracks.reserve(n_seeds);

  // Toal number of track states = number of tracks * max(states/track)
  const unsigned int n_max_states{n_seeds * cfg.max_track_candidates_per_track};

  // Create track states buffer (to make sure the resulting device container
  // is resizable)
  vecmem::copy copy{};
  typename edm::track_state_collection<algebra_t>::buffer track_states_buffer{
      cfg.run_smoother == smoother_type::e_none ? 0u : n_max_states, mr,
      vecmem::data::buffer_type::resizable};
  copy.setup(track_states_buffer)->ignore();

  // Track data collected by the measurement updater during pattern recog.
  vecmem::vector<track_state_candidate> track_cands{};
  vecmem::vector<filtered_track_state_candidate<algebra_t>>
      filtered_track_cands{};
  vecmem::vector<full_track_state_candidate<algebra_t>> full_track_cands{};

  switch (cfg.run_smoother) {
    case smoother_type::e_none: {
      track_cands.resize(n_max_states);
      break;
    }
    case smoother_type::e_kalman: {
      filtered_track_cands.resize(n_max_states);
      break;
    }
    case smoother_type::e_mbf: {
      full_track_cands.resize(n_max_states);
      break;
    }
    default: {
      TRACCC_FATAL_HOST("Unknown smoother option");
      return track_container;
    }
  }

  // Setup the surface sequence buffer
  const unsigned int n_surfaces_per_track{
      std::max(cfg.max_track_candidates_per_track *
                   cfg.kalman_smoother.surface_sequence_size_factor,
               cfg.kalman_smoother.min_surface_sequence_capacity)};
  std::vector<unsigned int> seqs_sizes(
      n_seeds,
      cfg.run_smoother == smoother_type::e_kalman ? n_surfaces_per_track : 0u);

  auto sf_sequences_buffer =
      vecmem::data::jagged_vector_buffer<typename detector_t::surface_type>{
          seqs_sizes, mr, &mr, vecmem::data::buffer_type::resizable};
  copy.setup(sf_sequences_buffer)->ignore();

  for (unsigned int seed_idx = 0u; seed_idx < seeds.size(); ++seed_idx) {
    const auto& seed = seeds[seed_idx];

    // Set the data pointer to the beginning of the range of the track
    const auto cand_offset{static_cast<unsigned int>(
        seed_idx * cfg.max_track_candidates_per_track)};

    track_state_candidate_data<algebra_t> candidate_data(
        cfg.run_smoother, cand_offset, vecmem::get_data(track_cands),
        vecmem::get_data(filtered_track_cands),
        vecmem::get_data(full_track_cands));

    // Run the progressive filter for this seed
    const track_stats<scalar_t> trk_stats =
        traccc::details::progressive_kalman_filter(
            det, field, measurements_view, vecmem::get_data(meas_ranges), seed,
            seed_idx, candidate_data.ptr(), sf_sequences_buffer, cfg);

    // Check track stats and build the new track object
    const unsigned int n_track_states{trk_stats.n_track_states};
    const int ndf_sum{static_cast<int>(trk_stats.ndf_sum) - 5};

    track_container.tracks.push_back({});

    // Bad track, don't write to container
    if (n_track_states == 0u) {
      continue;
    }

    edm::track track =
        track_container.tracks.at(track_container.tracks.size() - 1u);

    track.fit_outcome() = cfg.run_smoother == smoother_type::e_none
                              ? track_fit_outcome::UNKNOWN
                              : track_fit_outcome::SUCCESS;

    track.params() = seed;
    track.ndf() = static_cast<scalar_t>(ndf_sum);
    track.chi2() = trk_stats.chi2_sum;
    track.pval() = prob(trk_stats.chi2_sum, static_cast<scalar_t>(ndf_sum));
    track.nholes() = static_cast<unsigned int>(trk_stats.n_holes);
    track.constituent_links().resize(n_track_states);

    const auto track_state_offset{
        static_cast<unsigned int>(track_container.states.size())};

    // Generate the track states for this track
    for (unsigned int state_idx = track_state_offset;
         state_idx < track_state_offset + n_track_states; state_idx++) {
      const unsigned int link_idx{state_idx - track_state_offset};

      TRACCC_DEBUG_HOST("Adding track state (local idx "
                        << link_idx << ", global idx " << state_idx << ")");

      // Intermediate type required to build a view
      traccc::track_state_from_candidate<algebra_t>(
          candidate_data.ptr(), cfg.run_smoother, link_idx, measurements, track,
          track_states_buffer);
    }

    TRACCC_DEBUG_HOST("Added track " << track_container.tracks.size() - 1
                                     << " to track container: " << track);
  }

  copy(track_states_buffer, track_container.states)->wait();

  TRACCC_INFO_HOST("Track finding finished");

  // Run a Kalman filter in back propagation mode to smoothe the tracks
  if (cfg.run_smoother == smoother_type::e_kalman) {
    TRACCC_INFO_HOST("Run Kalman Smoother");

    using fitter_t = traccc::details::kalman_fitter_t<detector_t, bfield_t>;

    // Run the backwards track fitting (forwards is covered by PKF)
    traccc::host::kalman_smoother<fitter_t>(
        cfg.kalman_smoother, det, field, track_container, sf_sequences_buffer);

    TRACCC_INFO_HOST("Kalman Smoother: Finished");
  } else if (cfg.run_smoother == smoother_type::e_mbf) {
    TRACCC_FATAL_HOST("MBF not implemented for progressive Kalman Filter!");
  }

  return track_container;
}

}  // namespace traccc::host::details
