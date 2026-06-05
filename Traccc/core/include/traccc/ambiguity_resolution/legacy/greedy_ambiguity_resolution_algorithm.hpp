/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// Greedy ambiguity resolution adapted from ACTS code

namespace traccc::legacy {

/// Evicts tracks that seem to be duplicates or fakes. This algorithm takes a
/// greedy approach in the sense that it will remove the track which looks "most
/// duplicate/fake" first and continues the same process with the rest. That
/// process continues until the final state conditions are met.
///
/// The implementation works as follows:
///  1) Calculate shared hits per track.
///  2) If the maximum shared hits criteria is met, we are done.
///     This is the configurable amount of shared hits we are ok with
///     in our experiment.
///  3) Else, remove the track with the highest relative shared hits (i.e.
///     shared hits / hits).
///  4) Back to square 1.
class greedy_ambiguity_resolution_algorithm
    : public algorithm<edm::track_container<default_algebra>::host(
          const edm::track_container<default_algebra>::host&)>,
      public messaging {

    public:
    struct config_t {

        config_t() {}

        /// Maximum amount of shared hits per track. One (1) means "no shared
        /// hit allowed".
        std::uint32_t maximum_shared_hits = 1;

        /// Maximum number of iterations.
        std::uint32_t maximum_iterations = 1000000;

        /// Minimum number of measurement to form a track.
        std::size_t n_measurements_min = 3;
    };

    struct state_t {
        std::size_t number_of_tracks{};

        /// For this whole comment section, track_index refers to the index of a
        /// track in the initial input container.
        ///
        /// There is no (track_id) in this algorithm, only (track_index).

        /// Associates each track_index with the track's p-value
        std::vector<traccc::scalar> track_pval;

        /// Associates each track_index to the track's (measurement_id)s list
        std::vector<std::vector<std::size_t>> measurements_per_track;

        /// Associates each measurement_id to a set of (track_index)es sharing
        /// it
        std::unordered_map<std::size_t, std::set<std::size_t>>
            tracks_per_measurement;

        /// Associates each track_index to its number of shared measurements
        /// (among other tracks)
        std::vector<std::size_t> shared_measurements_per_track;

        /// Keeps the selected tracks indexes that have not (yet) been removed
        /// by the algorithm
        std::set<std::pair<std::size_t, std::size_t>> selected_tracks;
    };

    /// Constructor for the greedy ambiguity resolution algorithm
    ///
    /// @param cfg  Configuration object
    // greedy_ambiguity_resolution_algorithm(const config_type& cfg) :
    // _config(cfg) {}
    greedy_ambiguity_resolution_algorithm(
        const config_t cfg, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone())
        : messaging(std::move(logger)), _config{cfg}, m_mr{mr} {}

    /// Run the algorithm
    ///
    /// @param track_states the container of the fitted track parameters
    /// @return the container without ambiguous tracks
    output_type operator()(const edm::track_container<default_algebra>::host&
                               tracks) const override;

    /// Get configuration
    config_t& get_config() { return _config; }

    private:
    /// Computes the initial state for the input data. This function accumulates
    /// information that will later be used to accelerate the ambiguity
    /// resolution.
    ///
    /// @param tracks The input track container
    /// @param state An empty state object which is expected to be default
    /// constructed.
    void compute_initial_state(
        const edm::track_container<default_algebra>::host& tracks,
        state_t& state) const;

    /// Updates the state iteratively by evicting one track after the other
    /// until the final state conditions are met.
    ///
    /// @param state A state object that was previously filled by the
    /// initialization.
    void resolve(state_t& state) const;

    config_t _config;
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::legacy
