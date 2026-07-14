/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/fitting/status_codes.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace traccc::host::details {

/// Specialization for Triplet-based track fitting
///
///
/// @param[in] fitter           The triplet fitter object to use on the track
/// candidates
/// @param[in] track_container  Input track container
/// @param[in] mr               Memory resource to use for the output container
/// @param[in] copy             Copy object
///
/// @return A container of the fitted track states
///
template <typename algebra_t, typename fitter_t>
typename edm::track_container<algebra_t>::host triplet_fitting(
    fitter_t& fitter,
    const typename edm::track_container<algebra_t>::const_view& track_container,
    vecmem::memory_resource& mr, vecmem::copy& copy) {

    // Get the input collections(s).
    const edm::measurement_collection::const_device measurements{
        track_container.measurements};
    const typename edm::track_collection<algebra_t>::const_device
        track_candidates{track_container.tracks};

    // Create the output container.
    typename edm::track_container<algebra_t>::host result{
        mr, track_container.measurements};

    // Iterate over the tracks.
    for (unsigned int i = 0; i < track_candidates.size(); ++i) {

        // make empty track
        result.tracks.push_back({});
        auto fitted_track = result.tracks.at(result.tracks.size() - 1);

        // Input measurement indices to fitting
        vecmem::vector<unsigned int> measurement_idx;

        // Loop over links (measurements) in track candidate
        for (const edm::track_constituent_link& link :
             track_candidates.constituent_links().at(i)) {

            assert(link.type ==
                   edm::track_constituent_link::measurement);  // all links are
                                                               // measurements

            measurement_idx.push_back(link.index);
        }

        vecmem::data::vector_buffer<detray::geometry::identifier> seqs_buffer{};
        copy.setup(seqs_buffer)->wait();

        // Make triplets
        fitter.make_triplets(measurement_idx, measurements);

        // Run fitter
        auto first_state = fitter.fit(fitted_track, measurements);

        // Save track states in result
        unsigned state_counter = 0;
        for (const edm::track_constituent_link& link :
             track_candidates.constituent_links().at(i)) {

            assert(link.type ==
                   edm::track_constituent_link::measurement);  // all links are
                                                               // measurements

            // make link in the track
            fitted_track.constituent_links().push_back(
                {edm::track_constituent_link::track_state,
                 static_cast<unsigned int>(result.states.size())});

            // use fitted state
            if (state_counter == 0) {
                result.states.push_back(first_state);
            }

            else {
                // make a state from the measurement
                result.states.push_back(
                    edm::make_track_state<algebra_t>(measurements, link.index));
            }
            ++state_counter;
        }
    }

    // Return the result container.
    return result;
}

}  // namespace traccc::host::details
