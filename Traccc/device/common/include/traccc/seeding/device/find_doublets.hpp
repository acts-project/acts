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

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

namespace traccc::device {

/// Function finding all of the spacepoint doublets
///
/// Based on the information collected by @c traccc::device::count_doublets it
/// can fill collection with the specific doublet pairs that exist in the event.
///
/// @param[in] globalIndex       The index of the current thread
/// @param[in] config            Seedfinder configuration
/// @param[in] spacepoints       All spacepoints in the event
/// @param[in] sp_view           The spacepoint grid to find doublets on
/// @param[in] dc_view           Collection with the number of doublets to find
/// @param[out] mb_doublets_view Collection of middle-bottom doublets
/// @param[out] mt_doublets_view Collection of middle-top doublets
///
TRACCC_HOST_DEVICE
inline void find_doublets(
    global_index_t globalIndex, const seedfinder_config& config,
    const edm::spacepoint_collection::const_view& spacepoints,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const doublet_counter_collection_types::const_view& dc_view,
    device_doublet_collection_types::view mb_doublets_view,
    device_doublet_collection_types::view mt_doublets_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/find_doublets.ipp"
