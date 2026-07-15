/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/device/prefix_sum_element.hpp"
#include "traccc/edm/device/doublet_counter.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

namespace traccc::device {

/// Function used for calculating the number of spacepoint doublets
///
/// The count is necessary for allocating the appropriate amount of memory
/// for storing the information of the candidates in a next step.
///
/// @param[in] globalIndex   The index of the current thread
/// @param[in] config        Seedfinder configuration
/// @param[in] spacepoints   All the spacepoints of the event
/// @param[in] sp_view       The spacepoint grid to count doublets on
/// @param[in] sp_ps_view    Prefix sum for iterating over the spacepoint grid
/// @param[out] doublet_view Collection storing the number of doublets for each
/// spacepoint
/// @param[out] nMidBot      Total number of middle-bottom doublets
/// @param[out] nMidTop      Total number of middle-top doublets
///
TRACCC_HOST_DEVICE
inline void count_doublets(
    global_index_t globalIndex, const seedfinder_config& config,
    const edm::spacepoint_collection::const_view& spacepoints,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const vecmem::data::vector_view<const prefix_sum_element_t>& sp_ps_view,
    doublet_counter_collection_types::view doublet_view, unsigned int& nMidBot,
    unsigned int& nMidTop);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/count_doublets.ipp"
