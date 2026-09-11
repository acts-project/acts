/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/fitting/device/fit_payload.hpp"
#include "traccc/fitting/status_codes.hpp"

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_helpers.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Payload for the @c traccc::device::fit_prelude function
struct fit_prelude_payload {
  /// Input track parameter IDs
  vecmem::data::vector_view<const unsigned int> track_indices;
  /// Input tracks
  edm::track_container<default_algebra>::const_view input_tracks;
  /// Output tracks
  edm::track_container<default_algebra>::view output_tracks;
  /// Output track liveness
  vecmem::data::vector_view<unsigned int> track_liveness;

};  // struct fit_prelude_payload

/// Function to prepare the fitting payloads for the fitting algorithm
TRACCC_HOST_DEVICE inline void fit_prelude(const global_index_t globalIndex,
                                           const fit_prelude_payload& payload) {
  edm::track_collection<default_algebra>::device tracks(
      payload.output_tracks.tracks);

  if (globalIndex >= tracks.size()) {
    return;
  }

  edm::track_state_collection<default_algebra>::device track_states(
      payload.output_tracks.states);

  vecmem::device_vector<const unsigned int> track_indices(
      payload.track_indices);
  vecmem::device_vector<unsigned int> track_liveness(payload.track_liveness);

  const unsigned int param_id = track_indices.at(globalIndex);

  edm::track track = tracks.at(param_id);

  const edm::track_collection<default_algebra>::const_device track_candidates{
      payload.input_tracks.tracks};
  const edm::track track_candidate = track_candidates.at(param_id);
  const auto track_candidate_constituent_links =
      track_candidate.constituent_links();
  const edm::measurement_collection::const_device measurements{
      payload.input_tracks.measurements};
  for (const edm::track_constituent_link& link :
       track_candidate_constituent_links) {
    assert(link.type == edm::track_constituent_link::measurement);
    const unsigned int track_state_index = track_states.push_back(
        edm::make_track_state<default_algebra>(measurements, link.index));
    track.constituent_links().push_back(
        {edm::track_constituent_link::track_state, track_state_index});
  }

  // TODO: Set other stuff in the header?
  track.params() = track_candidate.params();
  track_liveness.at(param_id) = 1u;
}
}  // namespace traccc::device
