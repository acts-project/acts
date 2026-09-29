/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
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
#include "traccc/seeding/detail/lin_circle.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Function creating the linearised circle of one middle-top doublet
///
/// The linearised circles are computed once here, so that
/// @c traccc::device::count_triplets and @c traccc::device::find_triplets
/// do not need to recompute them.
///
/// @param[in] tid                 The index of the current thread
/// @param[in] mt_doublet_view     Collection storing the midTop doublets
/// @param[in] doublet_count_view  Collection of doublet counters
/// @param[in] spacepoint_view     All spacepoints in the event
/// @param[in] sp_grid_view        The spacepoint grid
/// @param[out] out_view           Linearised circles of the midTop doublets
///
TRACCC_HOST_DEVICE
inline void make_mid_top_lincircles(
    global_index_t tid,
    device::device_doublet_collection_types::const_view mt_doublet_view,
    device::doublet_counter_collection_types::const_view doublet_count_view,
    edm::spacepoint_collection::const_view spacepoint_view,
    traccc::details::spacepoint_grid_types::const_view sp_grid_view,
    vecmem::data::vector_view<lin_circle> out_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/make_mid_top_lincircles.ipp"
