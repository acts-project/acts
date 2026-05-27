/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/detail/seed_finding.hpp"

#include "doublet_finding.hpp"
#include "seed_filtering.hpp"
#include "triplet_finding.hpp"

namespace traccc::host::details {

struct seed_finding::impl {
    /// Constructor
    impl(const seedfinder_config& finder_config,
         const seedfilter_config& filter_config, vecmem::memory_resource& mr,
         std::unique_ptr<const Logger> logger)
        : m_midBot_finding(finder_config, mr,
                           logger->cloneWithSuffix("MidBotAlg")),
          m_midTop_finding(finder_config, mr,
                           logger->cloneWithSuffix("MidTopAlg")),
          m_triplet_finding(finder_config, filter_config, mr,
                            logger->cloneWithSuffix("TripletAlg")),
          m_seed_filtering(finder_config, filter_config, mr,
                           logger->cloneWithSuffix("FilterAlg")),
          m_mr{mr} {}

    /// Algorithm performing the mid bottom doublet finding
    doublet_finding<traccc::details::spacepoint_type::bottom> m_midBot_finding;
    /// Algorithm performing the mid top doublet finding
    doublet_finding<traccc::details::spacepoint_type::top> m_midTop_finding;
    /// Algorithm performing the triplet finding
    triplet_finding m_triplet_finding;
    /// Algorithm performing the seed selection
    seed_filtering m_seed_filtering;
    /// The memory resource to use
    vecmem::memory_resource& m_mr;
};

seed_finding::seed_finding(const seedfinder_config& finder_config,
                           const seedfilter_config& filter_config,
                           vecmem::memory_resource& mr,
                           std::unique_ptr<const Logger> logger)
    : messaging(logger->clone()),
      m_impl{std::make_unique<impl>(finder_config, filter_config, mr,
                                    std::move(logger))} {}

seed_finding::seed_finding(seed_finding&&) noexcept = default;

seed_finding::~seed_finding() = default;

seed_finding& seed_finding::operator=(seed_finding&&) noexcept = default;

edm::seed_collection::host seed_finding::operator()(
    const edm::spacepoint_collection::const_view& sp_view,
    const traccc::details::spacepoint_grid_types::host& sp_grid) const {

    // Create the result collection.
    edm::seed_collection::host seeds{m_impl->m_mr};

    // Create a device container for the spacepoints.
    const edm::spacepoint_collection::const_device spacepoints{sp_view};

    // Iterate over the spacepoint grid's bins.
    for (unsigned int i = 0; i < sp_grid.nbins(); ++i) {

        // Consider all spacepoints in this bin as "middle" spacepoints in the
        // seed.
        const auto& middle_indices = sp_grid.bin(i);

        // Evaluate these middle spacepoints one-by-one.
        for (unsigned int j = 0; j < middle_indices.size(); ++j) {

            // Internal identifier for this middle spacepoint.
            sp_location spM_location({i, j});

            // middule-bottom doublet search
            const auto mid_bot =
                m_impl->m_midBot_finding(spacepoints, sp_grid, spM_location);

            if (mid_bot.first.empty()) {
                continue;
            }

            // middule-top doublet search
            const auto mid_top =
                m_impl->m_midTop_finding(spacepoints, sp_grid, spM_location);

            if (mid_top.first.empty()) {
                continue;
            }

            triplet_collection_types::host triplets{&m_impl->m_mr};

            // triplet search from the combinations of two doublets which
            // share middle spacepoint
            for (unsigned int k = 0; k < mid_bot.first.size(); ++k) {

                const doublet& mid_bot_doublet = mid_bot.first[k];
                const lin_circle& mid_bot_lc = mid_bot.second[k];

                const triplet_collection_types::host triplets_for_mid_bot =
                    m_impl->m_triplet_finding(spacepoints, sp_grid,
                                              mid_bot_doublet, mid_bot_lc,
                                              mid_top.first, mid_top.second);

                triplets.insert(triplets.end(), triplets_for_mid_bot.begin(),
                                triplets_for_mid_bot.end());
            }

            // seed filtering
            m_impl->m_seed_filtering(spacepoints, sp_grid, triplets, seeds);
        }
    }

    return seeds;
}

}  // namespace traccc::host::details
