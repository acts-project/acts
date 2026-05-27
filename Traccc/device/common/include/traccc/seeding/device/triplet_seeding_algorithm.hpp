/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"
#include "traccc/device/prefix_sum_element.hpp"
#include "traccc/edm/device/device_doublet.hpp"
#include "traccc/edm/device/device_triplet.hpp"
#include "traccc/edm/device/doublet_counter.hpp"
#include "traccc/edm/device/triplet_counter.hpp"

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// System include(s).
#include <memory>

namespace traccc::device {

/// Main algorithm for performing triplet track seeding
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class triplet_seeding_algorithm
    : public algorithm<edm::seed_collection::buffer(
          const edm::spacepoint_collection::const_view&)>,
      public messaging,
      public algorithm_base {

    public:
    /// Constructor for the seed finding algorithm
    ///
    /// @param finder_config The seed finding configuration
    /// @param grid_config   The spacepoint grid forming configuration
    /// @param filter_config The seed filtering configuration
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param logger The logger instance to use
    ///
    triplet_seeding_algorithm(
        const seedfinder_config& finder_config,
        const spacepoint_grid_config& grid_config,
        const seedfilter_config& filter_config, const memory_resource& mr,
        vecmem::copy& copy,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());
    /// Destructor
    virtual ~triplet_seeding_algorithm();

    /// Operator executing the algorithm.
    ///
    /// @param spacepoints is a view of all spacepoints in the event
    /// @return the buffer of track seeds reconstructed from the spacepoints
    ///
    output_type operator()(const edm::spacepoint_collection::const_view&
                               spacepoints) const override;

    protected:
    /// @name Function(s) to be implemented by derived classes
    /// @{

    /// Payload for the @c count_grid_capacities_kernel function
    struct count_grid_capacities_kernel_payload {
        /// The number of spacepoints in the event
        edm::spacepoint_collection::const_view::size_type n_spacepoints;
        /// The seed finding configuration
        const seedfinder_config& config;
        /// The phi axis of the spacepoint grid
        const traccc::details::spacepoint_grid_types::host::axis_p0_type&
            phi_axis;
        /// The z axis of the spacepoint grid
        const traccc::details::spacepoint_grid_types::host::axis_p1_type&
            z_axis;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The buffer to write the grid capacities into
        vecmem::data::vector_view<unsigned int>& grid_capacities;
    };

    /// Spacepoint grid capacity counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void count_grid_capacities_kernel(
        const count_grid_capacities_kernel_payload& payload) const = 0;

    /// Payload for the @c populate_grid_kernel function
    struct populate_grid_kernel_payload {
        /// The number of spacepoints in the event
        edm::spacepoint_collection::const_view::size_type n_spacepoints;
        /// The seed finding configuration
        const seedfinder_config& config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The spacepoint grid to populate
        traccc::details::spacepoint_grid_types::view& grid;
        /// A prefix sum describing the grid contents
        vecmem::data::vector_view<prefix_sum_element_t>& grid_prefix_sum;
    };

    /// Spacepoint grid population kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void populate_grid_kernel(
        const populate_grid_kernel_payload& payload) const = 0;

    /// Payload for the @c count_doublets_kernel function
    struct count_doublets_kernel_payload {
        /// The number of spacepoints in the event
        edm::spacepoint_collection::const_view::size_type n_spacepoints;
        /// The seed finding configuration
        const seedfinder_config& config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The populated spacepoint grid
        const traccc::details::spacepoint_grid_types::const_view& grid;
        /// A prefix sum describing the grid contents
        const vecmem::data::vector_view<const prefix_sum_element_t>&
            grid_prefix_sum;
        /// The doublet counter collection to fill
        doublet_counter_collection_types::view& doublet_counter;
        /// The number of middle-bottom doublets found
        unsigned int& nMidBot;
        /// The number of middle-top doublets found
        unsigned int& nMidTop;
    };

    /// Doublet counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void count_doublets_kernel(
        const count_doublets_kernel_payload& payload) const = 0;

    /// Payload for the @c find_doublets_kernel function
    struct find_doublets_kernel_payload {
        /// The number of doublets counted earlier
        device::doublet_counter_collection_types::const_view::size_type
            n_doublets;
        /// The seed finding configuration
        const seedfinder_config& config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The populated spacepoint grid
        const traccc::details::spacepoint_grid_types::const_view& grid;
        /// The doublet counter collection
        const doublet_counter_collection_types::const_view& doublet_counter;
        /// The middle-bottom doublet collection to fill
        device_doublet_collection_types::view& mb_doublets;
        /// The middle-top doublet collection to fill
        device_doublet_collection_types::view& mt_doublets;
    };

    /// Doublet finding kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void find_doublets_kernel(
        const find_doublets_kernel_payload& payload) const = 0;

    /// Payload for the @c count_triplets_kernel function
    struct count_triplets_kernel_payload {
        /// The number of middle-bottom doublets found earlier
        unsigned int nMidBot;
        /// The seed finding configuration
        const seedfinder_config& config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The populated spacepoint grid
        const traccc::details::spacepoint_grid_types::const_view& grid;
        /// The doublet counter collection
        const doublet_counter_collection_types::const_view& doublet_counter;
        /// The middle-bottom doublet collection
        const device_doublet_collection_types::const_view& mb_doublets;
        /// The middle-top doublet collection
        const device_doublet_collection_types::const_view& mt_doublets;
        /// The triplet counter per middle spacepoint to fill
        triplet_counter_spM_collection_types::view& spM_counter;
        /// The triplet counter per middle-bottom doublet to fill
        triplet_counter_collection_types::view& midBot_counter;
    };

    /// Triplet counting kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void count_triplets_kernel(
        const count_triplets_kernel_payload& payload) const = 0;

    /// Payload for the @c triplet_counts_reduction_kernel function
    struct triplet_counts_reduction_kernel_payload {
        /// The number of doublets found earlier
        device::doublet_counter_collection_types::const_view::size_type
            n_doublets;
        /// The doublet counter collection
        const doublet_counter_collection_types::const_view& doublet_counter;
        /// The triplet counter per middle spacepoint
        triplet_counter_spM_collection_types::view& spM_counter;
        /// The total number of triplets found
        unsigned int& nTriplets;
    };

    /// Triplet count reduction kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void triplet_counts_reduction_kernel(
        const triplet_counts_reduction_kernel_payload& payload) const = 0;

    /// Payload for the @c find_triplets_kernel function
    struct find_triplets_kernel_payload {
        /// The number of middle-bottom doublets found earlier
        unsigned int nMidBot;
        /// The seed finding configuration
        const seedfinder_config& finding_config;
        /// The seed filtering configuration
        const seedfilter_config& filter_config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The populated spacepoint grid
        const traccc::details::spacepoint_grid_types::const_view& grid;
        /// The doublet counter collection
        const doublet_counter_collection_types::const_view& doublet_counter;
        /// The middle-top doublet collection
        const device_doublet_collection_types::const_view& mt_doublets;
        /// The triplet counter per middle spacepoint
        const triplet_counter_spM_collection_types::const_view& spM_tc;
        /// The triplet counter per middle-bottom doublet
        const triplet_counter_collection_types::const_view& midBot_tc;
        /// The triplet collection to fill
        device_triplet_collection_types::view& triplets;
    };

    /// Triplet finding kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void find_triplets_kernel(
        const find_triplets_kernel_payload& payload) const = 0;

    /// Payload for the @c update_triplet_weights_kernel function
    struct update_triplet_weights_kernel_payload {
        /// The number of triplets found earlier
        device_triplet_collection_types::const_view::size_type n_triplets;
        /// The seed filtering configuration
        const seedfilter_config& config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The triplet counter per middle spacepoint
        const triplet_counter_spM_collection_types::const_view& spM_tc;
        /// The triplet counter per middle-bottom doublet
        const triplet_counter_collection_types::const_view& midBot_tc;
        /// The triplet collection to update
        device_triplet_collection_types::view& triplets;
    };

    /// Triplet weight updater/filler kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void update_triplet_weights_kernel(
        const update_triplet_weights_kernel_payload& payload) const = 0;

    /// Payload for the @c select_seeds_kernel function
    struct select_seeds_kernel_payload {
        /// The number of doublets found earlier
        device::doublet_counter_collection_types::const_view::size_type
            n_doublets;
        /// The seed finding configuration
        const seedfinder_config& finder_config;
        /// The seed filtering configuration
        const seedfilter_config& filter_config;
        /// All spacepoints in the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The populated spacepoint grid
        const traccc::details::spacepoint_grid_types::const_view& grid;
        /// The triplet counter per middle spacepoint
        const triplet_counter_spM_collection_types::const_view& spM_tc;
        /// The triplet counter per middle-bottom doublet
        const triplet_counter_collection_types::const_view& midBot_tc;
        /// The triplet collection
        const device_triplet_collection_types::const_view& triplets;
        /// The seed collection to fill
        edm::seed_collection::view& seeds;
    };

    /// Seed selection/filling kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void select_seeds_kernel(
        const select_seeds_kernel_payload& payload) const = 0;

    /// @}

    private:
    /// Internal data type
    struct data;
    /// Pointer to internal data
    std::unique_ptr<data> m_data;

};  // class triplet_seeding_algorithm

}  // namespace traccc::device
