/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seed_finding.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_binning.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <memory>

namespace traccc::host {

/// Main algorithm for performing the track seeding on the CPU
class seeding_algorithm : public algorithm<edm::seed_collection::host(
                              const edm::spacepoint_collection::const_view&)>,
                          public messaging {

    public:
    /// Constructor for the seed finding algorithm
    ///
    /// @param finder_config The configuration for the seed finder
    /// @param grid_config The configuration for the spacepoint grid
    /// @param filter_config The configuration for the seed filter
    /// @param mr The memory resource to use
    ///
    seeding_algorithm(
        const seedfinder_config& finder_config,
        const spacepoint_grid_config& grid_config,
        const seedfilter_config& filter_config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Operator executing the algorithm.
    ///
    /// @param spacepoints All spacepoints in the event
    /// @return The track seeds reconstructed from the spacepoints
    ///
    output_type operator()(const edm::spacepoint_collection::const_view&
                               spacepoints) const override;

    private:
    /// Tool performing the spacepoint binning
    details::spacepoint_binning m_binning;
    /// Tool performing the seed finding
    details::seed_finding m_finding;

};  // class seeding_algorithm

}  // namespace traccc::host
