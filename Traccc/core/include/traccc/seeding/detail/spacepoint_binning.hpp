/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/utils/messaging.hpp"

// System include(s).
#include <functional>
#include <utility>

namespace traccc::host::details {

/// Spacepoint Binning for the seeding algorithm
class spacepoint_binning : public messaging {

    public:
    /// Constructor for the spacepoint binning
    ///
    /// @param config is seed finder configuration parameters
    /// @param grid_config is for spacepoint grid parameter
    /// @param mr is the vecmem memory resource
    ///
    spacepoint_binning(
        const seedfinder_config& config,
        const spacepoint_grid_config& grid_config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Operator executing the algorithm
    ///
    /// @param spacepoints All of the spacepoints of the event
    /// @return The spacepoints arranged in a Phi-Z grid
    ///
    traccc::details::spacepoint_grid_types::host operator()(
        const edm::spacepoint_collection::const_view& spacepoints) const;

    private:
    /// @name Tool configuration
    /// @{
    seedfinder_config m_config;
    spacepoint_grid_config m_grid_config;
    std::pair<traccc::details::spacepoint_grid_types::host::axis_p0_type,
              traccc::details::spacepoint_grid_types::host::axis_p1_type>
        m_axes;
    /// @}

    /// Memory resource to use
    std::reference_wrapper<vecmem::memory_resource> m_mr;

};  // class spacepoint_binning

}  // namespace traccc::host::details
