/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Payload structure for the @c device::gather_measurement_votes function
struct gather_measurement_votes_payload {
    vecmem::data::vector_view<const unsigned long long int> insertion_mutex;
    vecmem::data::vector_view<const unsigned int> tip_index;
    vecmem::data::vector_view<unsigned int> votes_per_tip;
    unsigned int max_num_tracks_per_measurement;
};

TRACCC_HOST_DEVICE inline void gather_measurement_votes(
    global_index_t thread_id, const gather_measurement_votes_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/gather_measurement_votes.ipp"
