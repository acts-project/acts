/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/device_triplet.hpp"
#include "traccc/edm/device/triplet_counter.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

namespace traccc::device {

/// Function used for updating the triplets' weights
///
/// @param[in] globalIndex   The index of the current thread
/// @param[in] filter_config Seedfilter configuration
/// @param[in] spacepoints   All spacepoints in the event
/// @param[in] sp_view       The spacepoint grid
/// @param[in] spM_tc_view   Collection of triplet counts per spM
/// @param[in] tc_view       Collection of triplet counts per midBot doublet
/// @param[in] data Array for temporary storage of quality parameters for
/// comparison of triplets
/// @param[inout] triplet_view Collection of triplets
///
TRACCC_HOST_DEVICE
inline void update_triplet_weights(
    global_index_t globalIndex, const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_view& spacepoints,
    const triplet_counter_spM_collection_types::const_view& spM_tc_view,
    const triplet_counter_collection_types::const_view& tc_view, scalar* data,
    device_triplet_collection_types::view triplet_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/update_triplet_weights.ipp"
