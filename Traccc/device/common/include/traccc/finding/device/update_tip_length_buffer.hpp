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

/// Payload for the @c device::update_tip_length_buffer function
struct update_tip_length_buffer_payload {
    vecmem::data::vector_view<const unsigned int> old_tip_length;
    vecmem::data::vector_view<unsigned int> new_tip_length;
    vecmem::data::vector_view<const unsigned int> measurement_votes;
    vecmem::data::vector_view<unsigned int> tip_to_output_map;
    float min_measurement_voting_fraction;
};

TRACCC_HOST_DEVICE inline void update_tip_length_buffer(
    global_index_t thread_id, const update_tip_length_buffer_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/update_tip_length_buffer.ipp"
