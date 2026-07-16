/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"

// System include
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <vector>

namespace traccc::host {

auto greedy_ambiguity_resolution_algorithm::operator()(
    const edm::track_container<default_algebra>::const_view& track_view) const
    -> output_type {

    const edm::track_container<default_algebra>::const_device input_tracks(
        track_view);

    const std::size_t n_tracks = input_tracks.tracks.size();

    // Make the output container
    edm::track_container<default_algebra>::host output{m_mr.get()};
    output.measurements = track_view.measurements;

    if (n_tracks == 0) {
        return output;
    }

    // Make sure that max_shared_meas is largen than zero
    assert(m_config.max_shared_meas > 0u);

    // Accepted ids to iterate
    std::vector<unsigned int> accepted_ids(n_tracks);
    std::iota(accepted_ids.begin(), accepted_ids.end(), 0);

    // Make measurement ID, pval and n_measurement vector
    std::vector<std::vector<std::size_t>> meas_ids(n_tracks);
    std::vector<traccc::scalar> pvals(n_tracks);
    std::vector<std::size_t> n_meas(n_tracks, 0);

    for (unsigned int i = 0; i < n_tracks; i++) {
        // Fill the pval vectors
        pvals[i] = input_tracks.tracks.at(i).pval();

        const vecmem::device_vector<const edm::track_constituent_link>
            measurement_links = input_tracks.tracks.at(i).constituent_links();
        const unsigned int n_cands = measurement_links.size();

        if (n_cands < m_config.min_meas_per_track) {
            // Reject if the number of measurements is less than the cut
            const auto it =
                std::lower_bound(accepted_ids.begin(), accepted_ids.end(), i);
            assert(it != accepted_ids.end() && *it == i);
            accepted_ids.erase(it);
        } else {
            // Fill measurement ids and n_measurements
            meas_ids[i].reserve(n_cands);
            for (const auto [type, idx] : measurement_links) {
                assert(type == edm::track_constituent_link::measurement);
                meas_ids[i].push_back(
                    input_tracks.measurements.at(idx).identifier());
            }
            n_meas[i] = n_cands;
            assert(n_cands == meas_ids[i].size());
        }
    }

    // Get the sorted unique measurement vector
    std::unordered_set<std::size_t> unique_meas_set;
    for (const unsigned int& i : accepted_ids) {
        const vecmem::device_vector<const edm::track_constituent_link>
            measurement_links = input_tracks.tracks.at(i).constituent_links();
        for (const auto [type, idx] : measurement_links) {
            assert(type == edm::track_constituent_link::measurement);
            unique_meas_set.insert(
                input_tracks.measurements.at(idx).identifier());
        }
    }

    std::vector<std::size_t> unique_meas(unique_meas_set.begin(),
                                         unique_meas_set.end());
    std::sort(unique_meas.begin(), unique_meas.end());

    // Record the tracks per measurement
    std::vector<std::vector<unsigned int>> tracks_per_measurement;
    tracks_per_measurement.resize(unique_meas.size());

    for (const unsigned int i : accepted_ids) {

        std::unordered_set<std::size_t> deduplicated_ids(meas_ids[i].begin(),
                                                         meas_ids[i].end());

        for (const std::size_t meas_id : deduplicated_ids) {

            const auto it = std::lower_bound(unique_meas.begin(),
                                             unique_meas.end(), meas_id);
            assert(it != unique_meas.end());
            const std::size_t unique_meas_idx = static_cast<std::size_t>(
                std::distance(unique_meas.begin(), it));
            tracks_per_measurement[unique_meas_idx].push_back(i);
        }
    }

    // Count the number of shared measurements
    std::vector<unsigned int> n_shared(n_tracks, 0);
    for (const unsigned int i : accepted_ids) {
        for (const std::size_t meas_id : meas_ids[i]) {
            const auto it = std::lower_bound(unique_meas.begin(),
                                             unique_meas.end(), meas_id);
            assert(it != unique_meas.end());
            const std::size_t unique_meas_idx = static_cast<std::size_t>(
                std::distance(unique_meas.begin(), it));
            if (tracks_per_measurement[unique_meas_idx].size() > 1) {
                n_shared[i]++;
            }
        }
    }

    // Make relative number of shared measurement vector
    std::vector<traccc::scalar> rel_shared(n_tracks);
    for (const unsigned int i : accepted_ids) {
        rel_shared[i] = static_cast<traccc::scalar>(n_shared[i]) /
                        static_cast<traccc::scalar>(n_meas[i]);
    }

    // Sort the track id with rel_shared and pval to find the worst track fast
    std::vector<unsigned int> sorted_ids = accepted_ids;

    auto track_comparator = [&rel_shared, &pvals](unsigned int a,
                                                  unsigned int b) {
        if (rel_shared[a] != rel_shared[b]) {
            return rel_shared[a] < rel_shared[b];
        }
        return pvals[a] > pvals[b];
    };
    std::sort(sorted_ids.begin(), sorted_ids.end(), track_comparator);

    // Iterate over tracks
    for (unsigned int iter = 0; iter < m_config.max_iterations; iter++) {
        // Terminate if there are no tracks to iterate
        if (accepted_ids.empty()) {
            break;
        }

        unsigned int max_shared(0);
        for (const unsigned int i : accepted_ids) {
            if (n_shared[i] > max_shared) {
                max_shared = n_shared[i];
            }
        }

        // Terminate if the max shared measurements is less than the cut value
        if (max_shared < m_config.max_shared_meas) {
            break;
        }

        // The last element of sorted vector is the worst track
        const unsigned int worst_track = sorted_ids.back();

        // Remove the worst track from the accepted ids
        const auto it1 = std::lower_bound(accepted_ids.begin(),
                                          accepted_ids.end(), worst_track);
        assert(it1 != accepted_ids.end() && *it1 == worst_track);
        accepted_ids.erase(it1);

        // Pop the worst (rejected) id from the sorted ids
        sorted_ids.pop_back();

        std::unordered_set<std::size_t> seen;
        std::vector<std::size_t> meas_ids_to_remove;
        for (const std::size_t id : meas_ids[worst_track]) {
            if (seen.insert(id).second) {
                meas_ids_to_remove.push_back(id);
            }
        }

        for (const std::size_t id : meas_ids_to_remove) {
            const auto it =
                std::lower_bound(unique_meas.begin(), unique_meas.end(), id);
            assert(it != unique_meas.end());
            const std::size_t unique_meas_idx = static_cast<std::size_t>(
                std::distance(unique_meas.begin(), it));

            std::vector<unsigned int>& tracks =
                tracks_per_measurement[unique_meas_idx];

            if (tracks.empty()) {
                continue;
            }

            // Remove the worst (rejected) id from the tracks associated
            // with measurement
            const auto it2 =
                std::lower_bound(tracks.begin(), tracks.end(), worst_track);
            assert(it2 != tracks.end() && *it2 == worst_track);
            tracks.erase(it2);

            // If there is only one track associated with measurement, the
            // number of shared measurement can be reduced by one
            if (tracks.size() == 1) {
                const unsigned int tid = tracks[0];

                n_shared[tid] -= static_cast<unsigned int>(
                    std::count(meas_ids[tid].begin(), meas_ids[tid].end(), id));
                rel_shared[tid] = static_cast<traccc::scalar>(n_shared[tid]) /
                                  static_cast<traccc::scalar>(n_meas[tid]);
                // Reposition the track to the next of the worse track which is
                // firstly found during the reverse iteration
                const auto it3 =
                    std::find(sorted_ids.begin(), sorted_ids.end(), tid);
                assert(it3 != sorted_ids.end());
                const auto it4 = std::lower_bound(sorted_ids.begin(), it3, tid,
                                                  track_comparator);

                if (it3 != it4) {
                    sorted_ids.erase(it3);
                    sorted_ids.insert(it4, tid);
                }
            }
        }

        // Make sure that sorted_ids stays sorted
        assert(std::is_sorted(sorted_ids.begin(), sorted_ids.end(),
                              track_comparator));
    }

    // Fill the output container with accepted tracks
    output.tracks.reserve(accepted_ids.size());
    for (const unsigned int i : accepted_ids) {
        // Get the track candidate proxy.
        const edm::track tcand = input_tracks.tracks.at(i);
        // Add it to the output container. In such a complicated way, because
        // the input container is a "device" container, and the output is a
        // "host" one. So pushing back the device proxy doesn't work directly.
        // (Jagged vectors cannot be converted to each other just so easily.)
        output.tracks.push_back({tcand.fit_outcome(),
                                 tcand.params(),
                                 tcand.ndf(),
                                 tcand.chi2(),
                                 tcand.pval(),
                                 tcand.nholes(),
                                 {tcand.constituent_links().begin(),
                                  tcand.constituent_links().end()}});
    }

    return output;
}

}  // namespace traccc::host
