/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/triplet_finding_helper.hpp"

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void find_triplets(
    const global_index_t globalIndex, const seedfinder_config& config,
    const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const doublet_counter_collection_types::const_view& dc_view,
    const device_doublet_collection_types::const_view& mid_top_doublet_view,
    const triplet_counter_spM_collection_types::const_view& spM_tc_view,
    const triplet_counter_collection_types::const_view& tc_view,
    device_triplet_collection_types::view triplet_view) {

    // Check if anything needs to be done.
    const triplet_counter_collection_types::const_device triplet_counts(
        tc_view);
    if (globalIndex >= triplet_counts.size()) {
        return;
    }

    // Get device copy of input parameters
    const edm::spacepoint_collection::const_device spacepoints{
        spacepoints_view};
    const doublet_counter_collection_types::const_device doublet_counts(
        dc_view);
    const device_doublet_collection_types::const_device mid_top_doublet_device(
        mid_top_doublet_view);
    const traccc::details::spacepoint_grid_types::const_device sp_grid(sp_view);
    const triplet_counter_spM_collection_types::const_device triplet_counts_spM(
        spM_tc_view);

    // Get the current work item information
    const triplet_counter mid_bot_counter = triplet_counts.at(globalIndex);
    const triplet_counter_spM spM_counter =
        triplet_counts_spM.at(mid_bot_counter.spM_counter_link);
    const doublet_counter doublet_count =
        doublet_counts.at(mid_bot_counter.spM_counter_link);

    const sp_location spM_loc = spM_counter.spM;
    const sp_location spB_loc = mid_bot_counter.spB;

    // middle spacepoint
    const unsigned int spM_idx = sp_grid.bin(spM_loc.bin_idx)[spM_loc.sp_idx];
    const edm::spacepoint_collection::const_device::const_proxy_type spM =
        spacepoints.at(spM_idx);

    // bottom spacepoint
    const unsigned int spB_idx = sp_grid.bin(spB_loc.bin_idx)[spB_loc.sp_idx];
    const edm::spacepoint_collection::const_device::const_proxy_type spB =
        spacepoints.at(spB_idx);

    // Set up the device result collection
    device_triplet_collection_types::device triplets(triplet_view);

    // Apply the conformal transformation to middle-bot doublet
    const traccc::lin_circle lb = doublet_finding_helper::transform_coordinates<
        details::spacepoint_type::bottom>(spM, spB);

    // Calculate some physical quantities required for triplet compatibility
    // check
    const scalar iSinTheta2 = 1 + lb.cotTheta() * lb.cotTheta();
    const scalar scatteringInRegion2 = config.maxScatteringAngle2 * iSinTheta2 *
                                       config.sigmaScattering *
                                       config.sigmaScattering;

    // These two quantities are used as output parameters in
    // triplet_finding_helper::isCompatible but their values are irrelevant
    scalar curvature, impact_parameter;

    // find the reference (start) index of the mid-top doublet collection
    // item vector, where the doublets are recorded
    const unsigned int mt_start_idx = doublet_count.m_posMidTop;
    const unsigned int mt_end_idx = mt_start_idx + doublet_count.m_nMidTop;
    // The position in which these triplets should be filled is the sum of the
    // position for all triplets which share the same middle spacepoint
    // and the one for those which also share the same bottom spacepoint.
    unsigned int posTriplets =
        mid_bot_counter.posTriplets + spM_counter.posTriplets;

    // iterate over mid-top doublets
    for (unsigned int i = mt_start_idx; i < mt_end_idx; ++i) {
        const sp_location spT_loc = mid_top_doublet_device[i].sp2;

        const unsigned int spT_idx =
            sp_grid.bin(spT_loc.bin_idx)[spT_loc.sp_idx];
        const edm::spacepoint_collection::const_device::const_proxy_type spT =
            spacepoints.at(spT_idx);

        // Apply the conformal transformation to middle-top doublet
        const traccc::lin_circle lt =
            doublet_finding_helper::transform_coordinates<
                details::spacepoint_type::top>(spM, spT);

        // Check if mid-bot and mid-top doublets can form a triplet
        if (triplet_finding_helper::isCompatible(
                spM, lb, lt, config, iSinTheta2, scatteringInRegion2, curvature,
                impact_parameter)) {

            // Add triplet to jagged vector
            triplets.at(posTriplets++) = device_triplet(
                {spB_idx, spM_idx, spT_idx, globalIndex, curvature,
                 -impact_parameter * filter_config.impactWeightFactor,
                 lb.Zo()});
        }
    }
}

}  // namespace traccc::device
