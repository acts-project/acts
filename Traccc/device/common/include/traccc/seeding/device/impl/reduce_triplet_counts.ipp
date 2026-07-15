/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void reduce_triplet_counts(
    const global_index_t globalIndex,
    const doublet_counter_collection_types::const_view& dc_view,
    triplet_counter_spM_collection_types::view spM_tc_view,
    unsigned int& num_triplets) {

    // Get device copy of input parameters
    triplet_counter_spM_collection_types::device spM_counts(spM_tc_view);
    // Check if anything needs to be done.
    if (globalIndex >= spM_counts.size()) {
        return;
    }

    const doublet_counter_collection_types::const_device doublet_counts(
        dc_view);

    // This should only ever be called with two collections of equal size
    assert(doublet_counts.size() == spM_counts.size());

    // Get triplet counter for this middle spacepoint
    triplet_counter_spM& this_spM_counter = spM_counts.at(globalIndex);

    // Check if anything needs to be done.
    if (this_spM_counter.m_nTriplets == 0) {
        return;
    }

    // Fill the middle spacepoint information of the spM triplet counter
    this_spM_counter.spM = doublet_counts.at(globalIndex).m_spM;

    // Increment total number of triplets and claim position for this middle
    // spacepoint's triplets
    vecmem::device_atomic_ref<unsigned int> nTriplets(num_triplets);
    this_spM_counter.posTriplets =
        nTriplets.fetch_add(this_spM_counter.m_nTriplets);
}

}  // namespace traccc::device
