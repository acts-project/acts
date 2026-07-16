/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/spacepoint_binning_helper.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void populate_grid(
    const global_index_t globalIndex, const seedfinder_config& config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    details::spacepoint_grid_types::view grid_view,
    vecmem::data::vector_view<prefix_sum_element_t> grid_prefix_sum_view) {

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

        // Set up the spacepoint grid object(s).
        details::spacepoint_grid_types::device grid(grid_view);
        const details::spacepoint_grid_types::device::axis_p0_type& phi_axis =
            grid.axis_p0();
        const details::spacepoint_grid_types::device::axis_p1_type& z_axis =
            grid.axis_p1();

        // Find the grid bin that the spacepoint belongs to.
        const unsigned int bin_index =
            phi_axis.bin(sp.phi()) + phi_axis.bins() * z_axis.bin(sp.z());

        // Add the spacepoint's index to the grid.
        const unsigned int sp_index =
            grid.bin(bin_index).push_back(globalIndex);

        // Add a prefix sum element for the spacepoint.
        vecmem::device_vector<prefix_sum_element_t> grid_prefix_sum(
            grid_prefix_sum_view);
        grid_prefix_sum.push_back({bin_index, sp_index});
    }
}

}  // namespace traccc::device
