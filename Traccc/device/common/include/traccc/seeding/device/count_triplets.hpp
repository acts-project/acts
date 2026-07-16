/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/device_doublet.hpp"
#include "traccc/edm/device/doublet_counter.hpp"
#include "traccc/edm/device/triplet_counter.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

namespace traccc::device {

/// Function used for calculating the number of spacepoint triplets
///
/// The count is necessary for allocating the appropriate amount of memory
/// for storing the information of the candidates in a next step.
///
/// @param[in] globalIndex          The index of the current thread
/// @param[in] config               Seedfinder configuration
/// @param[in] spacepoints          All spacepoints in the event
/// @param[in] sp_view              The spacepoint grid to count triplets on
/// @param[in] dc_view              Collection of doublet counters
/// @param[in] mid_bot_doublet_view Collection storing the midBot doublets
/// @param[in] mid_top_doublet_view Collection storing the midTop doublets
/// @param[out] spM_tc Collection storing the number of triplets per middle
/// spacepoint
/// @param[out] mb_tc  Collection storing the number of triplets per midBottom
/// doublet
///
TRACCC_HOST_DEVICE
inline void count_triplets(
    global_index_t globalIndex, const seedfinder_config& config,
    const edm::spacepoint_collection::const_view& spacepoints,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const doublet_counter_collection_types::const_view& dc_view,
    const device_doublet_collection_types::const_view& mid_bot_doublet_view,
    const device_doublet_collection_types::const_view& mid_top_doublet_view,
    triplet_counter_spM_collection_types::view spM_tc,
    triplet_counter_collection_types::view mb_tc);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/count_triplets.ipp"
