/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/spacepoint_binning_helper.hpp"

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void count_grid_capacities(
    const global_index_t globalIndex, const seedfinder_config& config,
    const details::spacepoint_grid_types::host::axis_p0_type& phi_axis,
    const details::spacepoint_grid_types::host::axis_p1_type& z_axis,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    vecmem::data::vector_view<unsigned int> grid_capacities_view) {

    // Check if anything needs to be done.
    const edm::spacepoint_collection::const_device spacepoints(
        spacepoints_view);
    if (globalIndex >= spacepoints.size()) {
        return;
    }
    const edm::spacepoint_collection::const_device::const_proxy_type sp =
        spacepoints.at(globalIndex);

    /// Check out if the spacepoint can be used for seeding.
    if (is_valid_sp(config, sp)) {

        // Find the grid bin that the spacepoint belongs to.
        const unsigned int bin_index =
            phi_axis.bin(sp.phi()) + phi_axis.bins() * z_axis.bin(sp.z());

        // Increase the capacity of the grid bin.
        vecmem::device_vector<unsigned int> grid_capacities(
            grid_capacities_view);
        vecmem::device_atomic_ref<unsigned int> bin_content(
            grid_capacities[bin_index]);
        bin_content.fetch_add(1);
    }
}

}  // namespace traccc::device
