/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "traccc/seeding/doublet_finding_helper.hpp"
#include "traccc/seeding/triplet_finding_helper.hpp"

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void count_triplets(
    const global_index_t globalIndex, const seedfinder_config& config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const doublet_counter_collection_types::const_view& dc_view,
    const device_doublet_collection_types::const_view& mid_bot_doublet_view,
    const device_doublet_collection_types::const_view& mid_top_doublet_view,
    triplet_counter_spM_collection_types::view spM_tc_view,
    triplet_counter_collection_types::view mb_tc_view) {

    // Create device copy of input parameters
    const device_doublet_collection_types::const_device mid_bot_doublet_device(
        mid_bot_doublet_view);
    // Check if anything needs to be done.
    if (globalIndex >= mid_bot_doublet_device.size()) {
        return;
    }

    // Get current mid bottom doublet
    const device_doublet mid_bot = mid_bot_doublet_device.at(globalIndex);

    // Create device copy of input parameters
    const edm::spacepoint_collection::const_device spacepoints{
        spacepoints_view};
    const device_doublet_collection_types::const_device mid_top_doublet_device(
        mid_top_doublet_view);
    const doublet_counter_collection_types::const_device dc_device(dc_view);

    // Create device copy of output parameterss
    triplet_counter_collection_types::device mb_triplet_counter(mb_tc_view);
    triplet_counter_spM_collection_types::device spM_triplet_counter(
        spM_tc_view);

    // Get all spacepoints
    const traccc::details::spacepoint_grid_types::const_device sp_device{
        sp_view};

    const unsigned int counter_link = mid_bot.counter_link;
    const doublet_counter doublet_counts = dc_device.at(counter_link);

    // middle spacepoint
    const sp_location spM_loc = doublet_counts.m_spM;
    const edm::spacepoint_collection::const_device::const_proxy_type spM =
        spacepoints.at(sp_device.bin(spM_loc.bin_idx)[spM_loc.sp_idx]);
    const sp_location spB_loc = mid_bot.sp2;
    // bottom spacepoint
    const edm::spacepoint_collection::const_device::const_proxy_type spB =
        spacepoints.at(sp_device.bin(spB_loc.bin_idx)[spB_loc.sp_idx]);

    // Apply the conformal transformation to middle-bot doublet
    traccc::lin_circle lb = doublet_finding_helper::transform_coordinates<
        details::spacepoint_type::bottom>(spM, spB);

    // Calculate some physical quantities required for triplet compatibility
    // check
    scalar iSinTheta2 = static_cast<scalar>(1.) + lb.cotTheta() * lb.cotTheta();
    scalar scatteringInRegion2 = config.maxScatteringAngle2 * iSinTheta2;
    scatteringInRegion2 *= config.sigmaScattering * config.sigmaScattering;

    // These two quantities are used as output parameters in
    // triplet_finding_helper::isCompatible but their values are irrelevant
    scalar curvature, impact_parameter;

    // find the reference (start) index of the mid-top doublet container
    // item vector, where the doublets are recorded
    const unsigned int mt_start_idx = doublet_counts.m_posMidTop;
    const unsigned int mt_end_idx = mt_start_idx + doublet_counts.m_nMidTop;

    // number of triplets per middle-bot doublet
    unsigned int num_triplets_per_mb = 0;

    // iterate over mid-top doublets
    for (unsigned int i = mt_start_idx; i < mt_end_idx; ++i) {
        const traccc::sp_location spT_loc = mid_top_doublet_device[i].sp2;

        const edm::spacepoint_collection::const_device::const_proxy_type spT =
            spacepoints.at(sp_device.bin(spT_loc.bin_idx)[spT_loc.sp_idx]);

        // Apply the conformal transformation to middle-top doublet
        traccc::lin_circle lt = doublet_finding_helper::transform_coordinates<
            details::spacepoint_type::top>(spM, spT);

        // Check if mid-bot and mid-top doublets can form a triplet
        if (triplet_finding_helper::isCompatible(
                spM, lb, lt, config, iSinTheta2, scatteringInRegion2, curvature,
                impact_parameter)) {
            num_triplets_per_mb++;
        }
    }

    // if the number of triplets per mb is larger than 0, write the triplet
    // counter into the collection
    if (num_triplets_per_mb > 0) {
        triplet_counter_spM& header = spM_triplet_counter.at(counter_link);
        vecmem::device_atomic_ref<unsigned int> nTriplets(header.m_nTriplets);
        const unsigned int posTriplets =
            nTriplets.fetch_add(num_triplets_per_mb);

        mb_triplet_counter.push_back(
            {spB_loc, counter_link, num_triplets_per_mb, posTriplets});
    }
}

}  // namespace traccc::device
