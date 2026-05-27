/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/fitting/status_codes.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

template <typename algebra_t>
TRACCC_HOST_DEVICE inline void fit_prelude(
    const global_index_t globalIndex,
    vecmem::data::vector_view<const unsigned int> param_ids_view,
    typename edm::track_container<algebra_t>::const_view track_candidates_view,
    typename edm::track_container<algebra_t>::view tracks_view,
    vecmem::data::vector_view<unsigned int> param_liveness_view) {

    typename edm::track_collection<algebra_t>::device tracks(
        tracks_view.tracks);

    if (globalIndex >= tracks.size()) {
        return;
    }

    typename edm::track_state_collection<algebra_t>::device track_states(
        tracks_view.states);

    vecmem::device_vector<const unsigned int> param_ids(param_ids_view);
    vecmem::device_vector<unsigned int> param_liveness(param_liveness_view);

    const unsigned int param_id = param_ids.at(globalIndex);

    edm::track track = tracks.at(param_id);

    const typename edm::track_collection<algebra_t>::const_device
        track_candidates{track_candidates_view.tracks};
    const edm::track track_candidate = track_candidates.at(param_id);
    const auto track_candidate_constituent_links =
        track_candidate.constituent_links();
    const edm::measurement_collection::const_device measurements{
        track_candidates_view.measurements};
    for (const edm::track_constituent_link& link :
         track_candidate_constituent_links) {
        assert(link.type == edm::track_constituent_link::measurement);
        const unsigned int track_state_index = track_states.push_back(
            edm::make_track_state<algebra_t>(measurements, link.index));
        track.constituent_links().push_back(
            {edm::track_constituent_link::track_state, track_state_index});
    }

    // TODO: Set other stuff in the header?
    track.params() = track_candidate.params();
    param_liveness.at(param_id) = 1u;
}

}  // namespace traccc::device
