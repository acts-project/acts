/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/finding/candidate_link.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Payload for the @c device::gather_best_tips_per_measurement function
template <typename algebra_t>
struct gather_best_tips_per_measurement_payload {
    vecmem::data::vector_view<const unsigned int> tips;
    vecmem::data::vector_view<const candidate_link> links;
    edm::measurement_collection::const_view measurements;
    vecmem::data::vector_view<unsigned long long int> insertion_mutex;
    vecmem::data::vector_view<unsigned int> tip_index;
    vecmem::data::vector_view<typename algebra_t::scalar> tip_pval;
    unsigned int max_num_tracks_per_measurement;
};

template <typename algebra_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gather_best_tips_per_measurement(
    global_index_t thread_id, const barrier_t& barrier,
    const gather_best_tips_per_measurement_payload<algebra_t>& payload);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/gather_best_tips_per_measurement.ipp"
