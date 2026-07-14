/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/ambiguity_resolution/legacy/greedy_ambiguity_resolution_algorithm.hpp"

// System include
#include <algorithm>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

// VecMem include(s).
#include <vecmem/containers/vector.hpp>

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/utils/algorithm.hpp"

// Greedy ambiguity resolution adapted from ACTS code

namespace traccc::legacy {

/// Run the algorithm
///
/// @param track_states the container of the fitted track parameters
/// @return the container without ambiguous tracks
auto greedy_ambiguity_resolution_algorithm::operator()(
    const edm::track_container<default_algebra>::host& track_states) const
    -> output_type {

    state_t state;
    compute_initial_state(track_states, state);
    resolve(state);

    // Copy the tracks to be retained in the return value

    edm::track_container<default_algebra>::host res{m_mr.get(),
                                                    track_states.measurements};
    res.tracks.reserve(state.selected_tracks.size());

    TRACCC_DEBUG(
        "state.selected_tracks.size() = " << state.selected_tracks.size());

    for (auto index : state.selected_tracks) {
        res.tracks.push_back(track_states.tracks.at(index.first));
    }
    return res;
}

void greedy_ambiguity_resolution_algorithm::compute_initial_state(
    const edm::track_container<default_algebra>::host& track_states,
    state_t& state) const {

    // Create a measurement collection to interact with.
    const edm::measurement_collection::const_device measurements(
        track_states.measurements);

    // For each track of the input container
    std::size_t n_track_states = track_states.tracks.size();
    for (std::size_t track_index = 0; track_index < n_track_states;
         ++track_index) {

        auto const track = track_states.tracks.at(track_index);

        // Kick out tracks that do not fulfill our initial requirements
        if (track.constituent_links().size() < _config.n_measurements_min) {
            continue;
        }

        // Create the list of measurement_id of the current track
        std::vector<std::size_t> measurement_ids;
        std::unordered_map<std::size_t, std::size_t> already_added_mes;

        for (const auto& [type, meas_idx] : track.constituent_links()) {
            assert(type == edm::track_constituent_link::measurement);

            std::size_t mid = measurements.at(meas_idx).identifier();

            // If the same measurement is found multiple times in a single
            // track: remove duplicates.
            auto dm = already_added_mes.insert(
                std::pair<std::size_t, std::size_t>(mid, 1));

            // If the measurement was already present in already_added_mes
            if (!dm.second) {
                // Increment the count for this measurement_id
                ++(dm.first->second);
                TRACCC_DEBUG("(1/3) Track " << track_index
                                            << " has duplicated measurement "
                                            << mid << ".");
            } else {
                measurement_ids.push_back(mid);
            }
        }

        // Add this track pval value
        state.track_pval.push_back(track.pval());
        // Add all the (measurement_id)s of this track
        state.measurements_per_track.push_back(std::move(measurement_ids));
        // Initially, every track is in the selected_track list. They will later
        // be removed according to the algorithm.
        state.selected_tracks.insert({track_index, state.number_of_tracks});
        ++state.number_of_tracks;
    }

    // Associate each measurement to the tracks sharing it
    for (std::size_t track_index = 0; track_index < state.number_of_tracks;
         ++track_index) {
        for (auto meas_id : state.measurements_per_track[track_index]) {
            state.tracks_per_measurement[meas_id].insert(track_index);
        }
    }

    // Finally, we can accumulate the number of shared measurements per track
    state.shared_measurements_per_track =
        std::vector<std::size_t>(state.number_of_tracks, 0);

    for (std::size_t track_index = 0; track_index < state.number_of_tracks;
         ++track_index) {
        for (auto meas_index : state.measurements_per_track[track_index]) {
            if (state.tracks_per_measurement[meas_index].size() > 1) {
                ++state.shared_measurements_per_track[track_index];
            }
        }
    }
}

namespace {

/// Removes a track from the state which has to be done for multiple properties
/// because of redundancy.
static void remove_track(greedy_ambiguity_resolution_algorithm::state_t& state,
                         std::pair<std::size_t, std::size_t> track_index) {
    for (auto meas_index : state.measurements_per_track[track_index.second]) {
        state.tracks_per_measurement[meas_index].erase(track_index.second);

        if (state.tracks_per_measurement[meas_index].size() == 1) {
            auto j_track = *state.tracks_per_measurement[meas_index].begin();
            --state.shared_measurements_per_track[j_track];
        }
    }
    state.selected_tracks.erase(track_index);
}
}  // namespace

void greedy_ambiguity_resolution_algorithm::resolve(state_t& state) const {
    /// Compares two tracks based on the number of shared measurements in order
    /// to decide if we already met the final state.
    auto shared_measurements_comperator =
        [&state](std::pair<std::size_t, std::size_t> a,
                 std::pair<std::size_t, std::size_t> b) {
            return state.shared_measurements_per_track[a.second] <
                   state.shared_measurements_per_track[b.second];
        };

    /// Compares two tracks in order to find the one which should be evicted.
    /// First we compare the relative amount of shared measurements. If that is
    /// indecisive we use the pval.
    auto track_comperator = [&state](std::pair<std::size_t, std::size_t> a,
                                     std::pair<std::size_t, std::size_t> b) {
        /// Helper to calculate the relative amount of shared measurements.
        auto relative_shared_measurements = [&state](std::size_t i) {
            return static_cast<double>(state.shared_measurements_per_track[i]) /
                   static_cast<double>(state.measurements_per_track[i].size());
        };

        if (relative_shared_measurements(a.second) !=
            relative_shared_measurements(b.second)) {
            return relative_shared_measurements(a.second) <
                   relative_shared_measurements(b.second);
        }
        return state.track_pval[a.second] > state.track_pval[b.second];
    };

    std::size_t iteration_count = 0;
    for (std::size_t i = 0; i < _config.maximum_iterations; ++i) {
        // Lazy out if there is nothing to filter on.
        if (state.selected_tracks.empty()) {
            TRACCC_DEBUG("No tracks left - exit loop");
            break;
        }

        // Find the maximum amount of shared measurements per track to decide if
        // we are done or not.
        auto maximum_shared_measurements = *std::max_element(
            state.selected_tracks.begin(), state.selected_tracks.end(),
            shared_measurements_comperator);

        if (state.shared_measurements_per_track[maximum_shared_measurements
                                                    .second] <
            _config.maximum_shared_hits) {
            break;
        }

        // Find the "worst" track by comparing them to each other
        auto bad_track =
            *std::max_element(state.selected_tracks.begin(),
                              state.selected_tracks.end(), track_comperator);

        remove_track(state, bad_track);
        ++iteration_count;
    }

    TRACCC_DEBUG("Iteration_count: " << iteration_count);
}

}  // namespace traccc::legacy
