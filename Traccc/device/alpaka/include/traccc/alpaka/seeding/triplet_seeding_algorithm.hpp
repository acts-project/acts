/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/seeding/device/triplet_seeding_algorithm.hpp"

namespace traccc::alpaka {

/// Main algorithm for performing the track seeding using Alpaka
class triplet_seeding_algorithm : public device::triplet_seeding_algorithm,
                                  public alpaka::algorithm_base {

    public:
    /// Constructor for the seed finding algorithm
    ///
    /// @param mr The memory resource to use
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param q  The Alpaka queue to perform the operations in
    ///
    triplet_seeding_algorithm(
        const seedfinder_config& finder_config,
        const spacepoint_grid_config& grid_config,
        const seedfilter_config& filter_config,
        const traccc::memory_resource& mr, vecmem::copy& copy, alpaka::queue& q,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    private:
    /// @name Function(s) inherited from @c traccc::device::seeding_algorithm
    /// @{

    /// Spacepoint grid capacity counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void count_grid_capacities_kernel(
        const count_grid_capacities_kernel_payload& payload) const override;

    /// Spacepoint grid population kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void populate_grid_kernel(
        const populate_grid_kernel_payload& payload) const override;

    /// Doublet counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void count_doublets_kernel(
        const count_doublets_kernel_payload& payload) const override;

    /// Doublet finding kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void find_doublets_kernel(
        const find_doublets_kernel_payload& payload) const override;

    /// Triplet counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void count_triplets_kernel(
        const count_triplets_kernel_payload& payload) const override;

    /// Triplet count reduction kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void triplet_counts_reduction_kernel(
        const triplet_counts_reduction_kernel_payload& payload) const override;

    /// Triplet finding kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void find_triplets_kernel(
        const find_triplets_kernel_payload& payload) const override;

    /// Triplet weight updater/filler kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void update_triplet_weights_kernel(
        const update_triplet_weights_kernel_payload& payload) const override;

    /// Seed selection/filling kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void select_seeds_kernel(
        const select_seeds_kernel_payload& payload) const override;

};  // class triplet_seeding_algorithm

}  // namespace traccc::alpaka
