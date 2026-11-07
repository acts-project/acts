/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/seed_collection.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Encoding of a seed weight that preserves its ordering as an unsigned
/// integer, with zero below every possible weight
TRACCC_HOST_DEVICE unsigned int encode_seed_weight(float weight);

/// Function recording, for every spacepoint of a seed, the highest weight of
/// the seeds it belongs to
///
/// @param[in] globalIndex        The index of the current thread
/// @param[in] seeds_view         Collection of seeds
/// @param[out] best_weight_view  Encoded best seed weight per spacepoint
///
TRACCC_HOST_DEVICE
inline void mark_seed_weights(
    global_index_t globalIndex,
    const edm::seed_collection::const_view& seeds_view,
    vecmem::data::vector_view<unsigned int> best_weight_view);

/// Function dropping the seeds whose three spacepoints all belong to seeds of
/// a higher weight
///
/// @param[in] globalIndex       The index of the current thread
/// @param[in] seeds_view        Collection of seeds
/// @param[in] best_weight_view  Encoded best seed weight per spacepoint
/// @param[out] out_view         Collection of the selected seeds
///
TRACCC_HOST_DEVICE
inline void suppress_seeds(
    global_index_t globalIndex,
    const edm::seed_collection::const_view& seeds_view,
    const vecmem::data::vector_view<const unsigned int>& best_weight_view,
    edm::seed_collection::view out_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/suppress_seeds.ipp"
