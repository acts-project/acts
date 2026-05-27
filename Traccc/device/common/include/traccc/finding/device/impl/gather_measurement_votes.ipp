/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/array_insertion_mutex.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE inline void gather_measurement_votes(
    global_index_t thread_id, const gather_measurement_votes_payload& payload) {

    const unsigned int measurement_idx =
        thread_id / payload.max_num_tracks_per_measurement;
    const unsigned int tip_idx =
        thread_id % payload.max_num_tracks_per_measurement;

    const vecmem::device_vector<const unsigned long long int> insertion_mutex(
        payload.insertion_mutex);
    const vecmem::device_vector<const unsigned int> tip_index(
        payload.tip_index);
    vecmem::device_vector<unsigned int> votes_per_tip(payload.votes_per_tip);

    if (measurement_idx >= insertion_mutex.size()) {
        return;
    }

    auto [locked, size, worst] =
        decode_insertion_mutex(insertion_mutex.at(measurement_idx));

    if (tip_idx < size) {
        vecmem::device_atomic_ref<unsigned int>(
            votes_per_tip.at(tip_index.at(thread_id)))
            .fetch_add(1u);
    }
}

}  // namespace traccc::device
