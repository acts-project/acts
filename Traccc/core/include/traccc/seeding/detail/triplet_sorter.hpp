/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/detail/triplet.hpp"

namespace traccc::details {

/// Functor to sort triplets with, during their final selection.
class triplet_sorter {

    public:
    /// Constructor
    ///
    /// @param[in] spacepoints All spacepoints in the event
    ///
    TRACCC_HOST_DEVICE
    triplet_sorter(const edm::spacepoint_collection::const_device& spacepoints,
                   const spacepoint_grid_types::const_device& sp_grid)
        : m_spacepoints{&spacepoints}, m_sp_grid{&sp_grid} {}

    /// Compare two triplets.
    ///
    /// @param[in] seed1 First triplet
    /// @param[in] seed2 Second triplet
    ///
    /// @return A boolean value to sort the triplets with
    ///
    TRACCC_HOST_DEVICE
    bool operator()(const triplet& seed1, const triplet& seed2) const {
        // If their weights are different, sort them just based on them.
        if (seed1.weight != seed2.weight) {
            return seed1.weight > seed2.weight;
        }
        // If they are the same, sort them based on the sum of the squares of
        // the y and z coordinates of the spacepoints.
        else {
            // Access the spacepoints.
            const edm::spacepoint_collection::const_device::const_proxy_type
                spB1 = m_spacepoints->at(
                    m_sp_grid->bin(seed1.sp1.bin_idx)[seed1.sp1.sp_idx]);
            const edm::spacepoint_collection::const_device::const_proxy_type
                spT1 = m_spacepoints->at(
                    m_sp_grid->bin(seed1.sp3.bin_idx)[seed1.sp3.sp_idx]);
            const edm::spacepoint_collection::const_device::const_proxy_type
                spB2 = m_spacepoints->at(
                    m_sp_grid->bin(seed2.sp1.bin_idx)[seed2.sp1.sp_idx]);
            const edm::spacepoint_collection::const_device::const_proxy_type
                spT2 = m_spacepoints->at(
                    m_sp_grid->bin(seed2.sp3.bin_idx)[seed2.sp3.sp_idx]);

            // Calculate these custom weights.
            const scalar seed1_sum = spB1.y() * spB1.y() + spB1.z() * spB1.z() +
                                     spT1.y() * spT1.y() + spT1.z() * spT1.z();
            const scalar seed2_sum = spB2.y() * spB2.y() + spB2.z() * spB2.z() +
                                     spT2.y() * spT2.y() + spT2.z() * spT2.z();

            // Compare the custom weights. These should be different in all
            // cases.
            return seed1_sum > seed2_sum;
        }
    }

    private:
    /// All spacepoints in the event
    const edm::spacepoint_collection::const_device* m_spacepoints;
    /// The spacepoint grid
    const spacepoint_grid_types::const_device* m_sp_grid;

};  // struct triplet_sorter

}  // namespace traccc::details
