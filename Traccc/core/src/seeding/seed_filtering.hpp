/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/detail/triplet.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>

namespace traccc::host::details {

/// Seed filtering to filter out the bad triplets
class seed_filtering : public messaging {

    public:
    /// Constructor with the seed filter configuration
    seed_filtering(
        const seedfinder_config& finding_config,
        const seedfilter_config& filter_config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Callable operator for the seed filtering
    ///
    /// @param[in] spacepoints All spacepoints in the event
    /// @param[in] sp_grid The spacepoint grid
    /// @param[in,out] triplets is the vector of triplets per middle spacepoint
    /// @param[out] seeds are the vector of seeds where the new compatible seeds
    ///             are added
    ///
    void operator()(const edm::spacepoint_collection::const_device& spacepoints,
                    const traccc::details::spacepoint_grid_types::host& sp_grid,
                    triplet_collection_types::host& triplets,
                    edm::seed_collection::host& seeds) const;

    private:
    /// Seed finder configuration
    seedfinder_config m_finder_config;
    /// Seed filter configuration
    seedfilter_config m_filter_config;
    /// The memory resource to use
    std::reference_wrapper<vecmem::memory_resource> m_mr;

};  // class seed_filtering

}  // namespace traccc::host::details
