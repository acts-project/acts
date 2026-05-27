/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <memory>

namespace traccc::host::details {

/// Tool for performing seed finding from binned spacepoints
class seed_finding : public messaging {

    public:
    /// Constructor for the seed finding
    ///
    /// @param find_config is seed finder configuration parameters
    /// @param filter_config is the seed filter configuration
    /// @param mr The memory resource to use
    ///
    seed_finding(
        const seedfinder_config& find_config,
        const seedfilter_config& filter_config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());
    /// Move constructor
    seed_finding(seed_finding&&) noexcept;
    /// Destructor
    ~seed_finding();

    /// Move assignment operator
    seed_finding& operator=(seed_finding&&) noexcept;

    /// Callable operator for the seed finding
    ///
    /// @param spacepoints All spacepoints in the event
    /// @param sp_grid The same spacepoints arranged in a 2D Phi-Z grid
    /// @return The spacepoint triplets that form the track seeds
    ///
    edm::seed_collection::host operator()(
        const edm::spacepoint_collection::const_view& spacepoints,
        const traccc::details::spacepoint_grid_types::host& sp_grid) const;

    private:
    /// Internal implementation struct
    struct impl;
    /// Pointer to the internal implementation
    std::unique_ptr<impl> m_impl;

};  // class seed_finding

}  // namespace traccc::host::details
