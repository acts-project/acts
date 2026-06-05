/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/detail/spacepoint_binning.hpp"

#include "traccc/seeding/spacepoint_binning_helper.hpp"

// Detray include(s).
#include <detray/definitions/indexing.hpp>

namespace traccc::host::details {

spacepoint_binning::spacepoint_binning(
    const seedfinder_config& config, const spacepoint_grid_config& grid_config,
    vecmem::memory_resource& mr, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_config(config),
      m_grid_config(grid_config),
      m_axes(get_axes(grid_config, mr)),
      m_mr(mr) {}

traccc::details::spacepoint_grid_types::host spacepoint_binning::operator()(
    const edm::spacepoint_collection::const_view& sp_view) const {

    // Set up a device container on top of the input.
    const edm::spacepoint_collection::const_device spacepoints{sp_view};

    // Create the result object.
    traccc::details::spacepoint_grid_types::host result{
        m_axes.first, m_axes.second, m_mr.get()};
    const auto& phi_axis = result.axis_p0();
    const auto& z_axis = result.axis_p1();

    // Arrange the spacepoints into the bins of the 2D grid.
    for (unsigned int i = 0; i < spacepoints.size(); ++i) {

        // Get a proxy for this spacepoint.
        const edm::spacepoint_collection::const_device::const_proxy_type sp =
            spacepoints.at(i);

        if (is_valid_sp(m_config, sp)) {
            const detray::dindex bin_index =
                phi_axis.bin(sp.phi()) + phi_axis.bins() * z_axis.bin(sp.z());
            result.bin(bin_index).push_back(i);
        }
    }
    return result;
}

}  // namespace traccc::host::details
