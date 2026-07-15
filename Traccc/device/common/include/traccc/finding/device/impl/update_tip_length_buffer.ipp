/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// System include(s).
#include <cassert>
#include <limits>

namespace traccc::device {

TRACCC_HOST_DEVICE inline void update_tip_length_buffer(
    global_index_t thread_id, const update_tip_length_buffer_payload& payload) {

    const vecmem::device_vector<const unsigned int> measurement_votes(
        payload.measurement_votes);
    if (thread_id >= measurement_votes.size()) {
        return;
    }

    const vecmem::device_vector<const unsigned int> old_tip_length(
        payload.old_tip_length);
    vecmem::device_vector<unsigned int> new_tip_length(payload.new_tip_length);

    const unsigned int total_measurements = old_tip_length.at(thread_id);
    const unsigned int total_votes = measurement_votes.at(thread_id);

    assert(total_votes <= total_measurements);

    const float vote_fraction = static_cast<float>(total_votes) /
                                static_cast<float>(total_measurements);

    vecmem::device_vector<unsigned int> tip_to_output_map(
        payload.tip_to_output_map);

    if (vote_fraction < payload.min_measurement_voting_fraction) {
        tip_to_output_map.at(thread_id) =
            std::numeric_limits<unsigned int>::max();
    } else {
        const vecmem::device_vector<unsigned int>::size_type new_idx =
            new_tip_length.push_back(total_measurements);
        tip_to_output_map.at(thread_id) = new_idx;
    }
}

}  // namespace traccc::device
