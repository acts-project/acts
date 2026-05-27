/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "seed_filtering.hpp"

#include "traccc/seeding/detail/triplet_sorter.hpp"
#include "traccc/seeding/seed_selecting_helper.hpp"

// System include(s).
#include <algorithm>

namespace traccc::host::details {

seed_filtering::seed_filtering(const seedfinder_config& finder_config,
                               const seedfilter_config& filter_config,
                               vecmem::memory_resource& mr,
                               std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      m_finder_config(finder_config),
      m_filter_config(filter_config),
      m_mr{mr} {}

void seed_filtering::operator()(
    const edm::spacepoint_collection::const_device& spacepoints,
    const traccc::details::spacepoint_grid_types::host& sp_grid,
    triplet_collection_types::host& triplets,
    edm::seed_collection::host& seeds) const {

    // Select the triplets passing the "single seed cuts".
    std::vector<std::reference_wrapper<const triplet>>
        triplets_passing_single_seed_cuts;
    triplets_passing_single_seed_cuts.reserve(triplets.size());
    for (triplet& triplet : triplets) {
        // bottom
        const sp_location& spB_location = triplet.sp1;
        const unsigned int spB_idx =
            sp_grid.bin(spB_location.bin_idx)[spB_location.sp_idx];
        const edm::spacepoint_collection::const_device::const_proxy_type spB =
            spacepoints.at(spB_idx);

        // middle
        const sp_location& spM_location = triplet.sp2;
        const unsigned int spM_idx =
            sp_grid.bin(spM_location.bin_idx)[spM_location.sp_idx];
        const edm::spacepoint_collection::const_device::const_proxy_type spM =
            spacepoints.at(spM_idx);

        // top
        const sp_location& spT_location = triplet.sp3;
        const unsigned int spT_idx =
            sp_grid.bin(spT_location.bin_idx)[spT_location.sp_idx];
        const edm::spacepoint_collection::const_device::const_proxy_type spT =
            spacepoints.at(spT_idx);

        // Updat the triplet weight in-situ.
        seed_selecting_helper::seed_weight(m_filter_config, spM, spB, spT,
                                           triplet.weight);

        // Check if the triplet passes the single seed cuts.
        if (!seed_selecting_helper::single_seed_cut(m_filter_config, spM, spB,
                                                    spT, triplet.weight)) {
            continue;
        }

        // If so, keep a pointer to it.
        triplets_passing_single_seed_cuts.push_back(triplet);
    }

    // sort seeds based on their weights
    const traccc::details::spacepoint_grid_types::const_data sp_grid_data =
        traccc::get_data(sp_grid, m_mr.get());
    std::sort(triplets_passing_single_seed_cuts.begin(),
              triplets_passing_single_seed_cuts.end(),
              traccc::details::triplet_sorter{spacepoints, sp_grid_data});

    // Select the best ones.
    std::vector<std::reference_wrapper<const triplet>>
        triplets_passing_final_cuts;
    triplets_passing_final_cuts.reserve(
        triplets_passing_single_seed_cuts.size());

    if (triplets_passing_single_seed_cuts.size() > 0u) {

        // Always accept the first element.
        triplets_passing_final_cuts.push_back(
            triplets_passing_single_seed_cuts[0]);

        // Consider only a maximum number of triplets for the final quality cut.
        const std::size_t itLength =
            std::min(triplets_passing_single_seed_cuts.size(),
                     static_cast<std::size_t>(m_finder_config.maxSeedsPerSpM));
        for (std::size_t i = 1; i < itLength; ++i) {
            const traccc::details::spacepoint_grid_types::const_device
                sp_grid_accessor(sp_grid_data);
            const auto& this_seed = triplets_passing_single_seed_cuts[i].get();
            if (seed_selecting_helper::cut_per_middle_sp(
                    m_filter_config,
                    spacepoints.at(sp_grid_accessor.bin(
                        this_seed.sp1.bin_idx)[this_seed.sp1.sp_idx]),
                    this_seed.weight)) {
                triplets_passing_final_cuts.push_back(
                    triplets_passing_single_seed_cuts[i]);
            }
        }
    }

    // Add the best remaining seeds to the output collection.
    for (std::size_t i = 0;
         const triplet& triplet : triplets_passing_final_cuts) {
        if (i++ >= m_finder_config.maxSeedsPerSpM) {
            break;
        }
        seeds.push_back({sp_grid.bin(triplet.sp1.bin_idx)[triplet.sp1.sp_idx],
                         sp_grid.bin(triplet.sp2.bin_idx)[triplet.sp2.sp_idx],
                         sp_grid.bin(triplet.sp3.bin_idx)[triplet.sp3.sp_idx],
                         static_cast<float>(triplet.weight)});
    }
}

}  // namespace traccc::host::details
