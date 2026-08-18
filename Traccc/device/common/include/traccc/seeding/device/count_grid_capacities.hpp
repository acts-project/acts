/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Function used for calculating the capacity for the spacepoint grid
///
/// Before filling the spacepoint grid with the spacepoints that belong to
/// each grid bin, we need to calculate how big each of those bins are going
/// to be.
///
/// This function needs to be called separately for every spacepoint of the
/// event.
///
/// @param[in] globalIndex   The index of the current thread
/// @param[in] config        Seedfinder configuration
/// @param[in] phi_axis      The circular &Phi axis describing the geometry
/// @param[in] z_axis        The linear Z axis describing the geometry
/// @param[in] spacepoints   All the spacepoints of the event
/// @param[out] grid_capacities Capacity required for each spacepoint grid bin
///
TRACCC_HOST_DEVICE
inline void count_grid_capacities(
    global_index_t globalIndex, const seedfinder_config& config,
    const details::spacepoint_grid_types::host::axis_p0_type& phi_axis,
    const details::spacepoint_grid_types::host::axis_p1_type& z_axis,
    const edm::spacepoint_collection::const_view& spacepoints,
    vecmem::data::vector_view<unsigned int> grid_capacities);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/count_grid_capacities.ipp"
