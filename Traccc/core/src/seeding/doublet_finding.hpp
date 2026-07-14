/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/seeding/detail/doublet.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/detail/spacepoint_type.hpp"
#include "traccc/seeding/doublet_finding_helper.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>
#include <utility>

namespace traccc::host::details {

/// Doublet finding to search the combinations of two compatible spacepoints
/// @tparam otherSpType is whether it is for middle-bottom or middle-top doublet
template <traccc::details::spacepoint_type otherSpType>
struct doublet_finding : public messaging {

    // A small sanity check
    static_assert(otherSpType == traccc::details::spacepoint_type::bottom ||
                  otherSpType == traccc::details::spacepoint_type::top);

    /// Constructor for the doublet finding
    ///
    /// @param config is the configuration parameters
    ///
    doublet_finding(
        const seedfinder_config& config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone())
        : messaging(std::move(logger)), m_config{config}, m_mr{mr} {}

    /// Callable operator for doublet finding per middle spacepoint
    ///
    /// @param spacepoints The spacepoint container
    /// @param sp_grid The spacepoint grid
    /// @param middle_sp The middle spacepoint to find doublets for
    /// @return A pair of vectors of doublets and transformed coordinates
    ///
    std::pair<doublet_collection_types::host, lin_circle_collection_types::host>
    operator()(const edm::spacepoint_collection::const_device& spacepoints,
               const traccc::details::spacepoint_grid_types::host& sp_grid,
               const sp_location& middle_location) const {

        // Create the result object.
        auto result =
            std::make_pair(doublet_collection_types::host{&(m_mr.get())},
                           lin_circle_collection_types::host{&(m_mr.get())});

        // Access the middle spacepoint.
        const edm::spacepoint_collection::const_device::const_proxy_type
            middle_sp = spacepoints.at(
                sp_grid.bin(middle_location.bin_idx)[middle_location.sp_idx]);

        // Get the Phi/Z bins in which to look for the other spacepoint of the
        // doublet.
        const detray::dindex_sequence phi_bins =
            sp_grid.axis_p0().zone(middle_sp.phi(), m_config.neighbor_scope);
        const detray::dindex_sequence z_bins =
            sp_grid.axis_p1().zone(middle_sp.z(), m_config.neighbor_scope);

        // Iterate over neighbor bins.
        for (detray::dindex phi_bin : phi_bins) {
            for (detray::dindex z_bin : z_bins) {

                // Get the global index for this bin.
                const detray::dindex bin_idx =
                    phi_bin + z_bin * sp_grid.axis_p0().bins();

                // Get the indices of the spacepoints in this bin.
                traccc::details::spacepoint_grid_types::host::
                    serialized_storage::const_reference sp_indices =
                        sp_grid.bin(phi_bin, z_bin);
                for (unsigned int i = 0; unsigned int sp_index : sp_indices) {

                    // Access the other spacepoint.
                    const edm::spacepoint_collection::const_device::
                        const_proxy_type other_sp = spacepoints.at(sp_index);

                    // Check if the spacepoints are compatible.
                    if (doublet_finding_helper::isCompatible<otherSpType>(
                            middle_sp, other_sp, m_config)) {

                        // If so, create a doublet for them.
                        result.first.push_back(
                            {middle_location,
                             {static_cast<unsigned int>(bin_idx), i}});
                        result.second.push_back(
                            doublet_finding_helper::transform_coordinates<
                                otherSpType>(middle_sp, other_sp));
                    }
                    ++i;
                }
            }
        }

        // Return the result.
        return result;
    }

    private:
    /// The doublet finding configuration parameters
    seedfinder_config m_config;
    /// The memory resource
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::host::details
