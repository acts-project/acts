/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/alpaka/seeding/triplet_seeding_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/seeding/device/count_doublets.hpp"
#include "traccc/seeding/device/count_grid_capacities.hpp"
#include "traccc/seeding/device/count_triplets.hpp"
#include "traccc/seeding/device/find_doublets.hpp"
#include "traccc/seeding/device/find_triplets.hpp"
#include "traccc/seeding/device/populate_grid.hpp"
#include "traccc/seeding/device/reduce_triplet_counts.hpp"
#include "traccc/seeding/device/select_seeds.hpp"
#include "traccc/seeding/device/update_triplet_weights.hpp"

namespace traccc::alpaka {
namespace kernels {

/// Kernel for running @c traccc::device::count_grid_capacities
struct count_grid_capacity {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const seedfinder_config& config,
        const traccc::details::spacepoint_grid_types::host::axis_p0_type
            phi_axis,
        const traccc::details::spacepoint_grid_types::host::axis_p1_type z_axis,
        const edm::spacepoint_collection::const_view spacepoints_view,
        vecmem::data::vector_view<unsigned int> grid_capacities_view) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::count_grid_capacities(globalThreadIdx, config, phi_axis, z_axis,
                                      spacepoints_view, grid_capacities_view);
    }
};

/// Kernel for running @c traccc::device::populate_grid
struct populate_grid {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config config,
        edm::spacepoint_collection::const_view spacepoints_view,
        traccc::details::spacepoint_grid_types::view grid_view,
        vecmem::data::vector_view<device::prefix_sum_element_t> grid_prefix_sum)
        const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::populate_grid(globalThreadIdx, config, spacepoints_view,
                              grid_view, grid_prefix_sum);
    }
};

/// Kernel for running @c traccc::device::count_doublets
struct count_doublets {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config config,
        const edm::spacepoint_collection::const_view spacepoints,
        const traccc::details::spacepoint_grid_types::const_view sp_grid,
        const vecmem::data::vector_view<const device::prefix_sum_element_t>
            sp_prefix_sum,
        device::doublet_counter_collection_types::view doublet_counter,
        unsigned int* nMidBot, unsigned int* nMidTop) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::count_doublets(globalThreadIdx, config, spacepoints, sp_grid,
                               sp_prefix_sum, doublet_counter, *nMidBot,
                               *nMidTop);
    }
};

/// Kernel for running @c traccc::device::find_doublets
struct find_doublets {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config config,
        edm::spacepoint_collection::const_view spacepoints,
        traccc::details::spacepoint_grid_types::const_view sp_grid,
        device::doublet_counter_collection_types::const_view doublet_counter,
        device::device_doublet_collection_types::view mb_doublets,
        device::device_doublet_collection_types::view mt_doublets) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::find_doublets(globalThreadIdx, config, spacepoints, sp_grid,
                              doublet_counter, mb_doublets, mt_doublets);
    }
};

/// Kernel for running @c traccc::device::count_triplets
struct count_triplets {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config config,
        edm::spacepoint_collection::const_view spacepoints,
        traccc::details::spacepoint_grid_types::const_view sp_grid,
        device::doublet_counter_collection_types::const_view doublet_counter,
        device::device_doublet_collection_types::const_view mb_doublets,
        device::device_doublet_collection_types::const_view mt_doublets,
        device::triplet_counter_spM_collection_types::view spM_counter,
        device::triplet_counter_collection_types::view midBot_counter) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::count_triplets(globalThreadIdx, config, spacepoints, sp_grid,
                               doublet_counter, mb_doublets, mt_doublets,
                               spM_counter, midBot_counter);
    }
};

/// Kernel for running @c traccc::device::reduce_triplet_counts
struct reduce_triplet_counts {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        device::doublet_counter_collection_types::const_view doublet_counter,
        device::triplet_counter_spM_collection_types::view spM_counter,
        unsigned int* nTriplets) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::reduce_triplet_counts(globalThreadIdx, doublet_counter,
                                      spM_counter, *nTriplets);
    }
};

/// Kernel for running @c traccc::device::find_triplets
struct find_triplets {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config config,
        seedfilter_config filter_config,
        edm::spacepoint_collection::const_view spacepoints,
        traccc::details::spacepoint_grid_types::const_view sp_grid,
        device::doublet_counter_collection_types::const_view doublet_counter,
        device::device_doublet_collection_types::const_view mt_doublets,
        device::triplet_counter_spM_collection_types::const_view spM_tc,
        device::triplet_counter_collection_types::const_view midBot_tc,
        device::device_triplet_collection_types::view triplet_view) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        device::find_triplets(globalThreadIdx, config, filter_config,
                              spacepoints, sp_grid, doublet_counter,
                              mt_doublets, spM_tc, midBot_tc, triplet_view);
    }
};

/// Kernel for running @c traccc::device::update_triplet_weights
struct update_triplet_weights {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfilter_config filter_config,
        edm::spacepoint_collection::const_view spacepoints,
        device::triplet_counter_spM_collection_types::const_view spM_tc,
        device::triplet_counter_collection_types::const_view midBot_tc,
        device::device_triplet_collection_types::view triplet_view) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        auto const localThreadIdx =
            ::alpaka::getIdx<::alpaka::Block, ::alpaka::Threads>(acc)[0u];

        // Array for temporary storage of quality parameters for comparing
        // triplets within weight updating kernel
        scalar* const data = ::alpaka::getDynSharedMem<scalar>(acc);

        // Each thread uses compatSeedLimit elements of the array
        scalar* dataPos = &data[localThreadIdx * filter_config.compatSeedLimit];

        device::update_triplet_weights(globalThreadIdx, filter_config,
                                       spacepoints, spM_tc, midBot_tc, dataPos,
                                       triplet_view);
    }
};

/// Kernel for running @c traccc::device::select_seeds
struct select_seeds {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, seedfinder_config finder_config,
        seedfilter_config filter_config,
        edm::spacepoint_collection::const_view spacepoints,
        traccc::details::spacepoint_grid_types::const_view sp_view,
        device::triplet_counter_spM_collection_types::const_view spM_tc,
        device::triplet_counter_collection_types::const_view midBot_tc,
        device::device_triplet_collection_types::const_view triplet_view,
        edm::seed_collection::view seed_view) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];
        auto const localThreadIdx =
            ::alpaka::getIdx<::alpaka::Block, ::alpaka::Threads>(acc)[0u];

        // Array for temporary storage of quality parameters for comparing
        // triplets within weight updating kernel
        device::device_triplet* const data =
            ::alpaka::getDynSharedMem<device::device_triplet>(acc);

        // Each thread uses max_triplets_per_spM elements of the array
        device::device_triplet* dataPos =
            &data[localThreadIdx * finder_config.maxSeedsPerSpM];

        device::select_seeds(globalThreadIdx, finder_config, filter_config,
                             spacepoints, sp_view, spM_tc, midBot_tc,
                             triplet_view, dataPos, seed_view);
    }
};

}  // namespace kernels

triplet_seeding_algorithm::triplet_seeding_algorithm(
    const seedfinder_config& finder_config,
    const spacepoint_grid_config& grid_config,
    const seedfilter_config& filter_config, const traccc::memory_resource& mr,
    vecmem::copy& copy, alpaka::queue& q, std::unique_ptr<const Logger> logger)
    : device::triplet_seeding_algorithm(finder_config, grid_config,
                                        filter_config, mr, copy,
                                        std::move(logger)),
      alpaka::algorithm_base{q} {}

void triplet_seeding_algorithm::count_grid_capacities_kernel(
    const count_grid_capacities_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::count_grid_capacity{}, payload.config, payload.phi_axis,
        payload.z_axis, payload.spacepoints, payload.grid_capacities);
}

void triplet_seeding_algorithm::populate_grid_kernel(
    const populate_grid_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::populate_grid{}, payload.config, payload.spacepoints,
        payload.grid, payload.grid_prefix_sum);
}

void triplet_seeding_algorithm::count_doublets_kernel(
    const count_doublets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::count_doublets{}, payload.config, payload.spacepoints,
        payload.grid, payload.grid_prefix_sum, payload.doublet_counter,
        &(payload.nMidBot), &(payload.nMidTop));
}

void triplet_seeding_algorithm::find_doublets_kernel(
    const find_doublets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::find_doublets{}, payload.config, payload.spacepoints,
        payload.grid, payload.doublet_counter, payload.mb_doublets,
        payload.mt_doublets);
}

void triplet_seeding_algorithm::count_triplets_kernel(
    const count_triplets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks = (payload.nMidBot + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::count_triplets{}, payload.config, payload.spacepoints,
        payload.grid, payload.doublet_counter, payload.mb_doublets,
        payload.mt_doublets, payload.spM_counter, payload.midBot_counter);
}

void triplet_seeding_algorithm::triplet_counts_reduction_kernel(
    const triplet_counts_reduction_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::reduce_triplet_counts{}, payload.doublet_counter,
        payload.spM_counter, &(payload.nTriplets));
}

void triplet_seeding_algorithm::find_triplets_kernel(
    const find_triplets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks = (payload.nMidBot + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::find_triplets{}, payload.finding_config, payload.filter_config,
        payload.spacepoints, payload.grid, payload.doublet_counter,
        payload.mt_doublets, payload.spM_tc, payload.midBot_tc,
        payload.triplets);
}

void triplet_seeding_algorithm::update_triplet_weights_kernel(
    const update_triplet_weights_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_triplets + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::update_triplet_weights{}, payload.config, payload.spacepoints,
        payload.spM_tc, payload.midBot_tc, payload.triplets);
}

void triplet_seeding_algorithm::select_seeds_kernel(
    const select_seeds_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), makeWorkDiv<Acc>(n_blocks, n_threads),
        kernels::select_seeds{}, payload.finder_config, payload.filter_config,
        payload.spacepoints, payload.grid, payload.spM_tc, payload.midBot_tc,
        payload.triplets, payload.seeds);
}

}  // namespace traccc::alpaka

// Define the required trait needed for Dynamic shared memory allocation.
namespace alpaka::trait {

template <typename TAcc>
struct BlockSharedMemDynSizeBytes<
    traccc::alpaka::kernels::update_triplet_weights, TAcc> {
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
        traccc::alpaka::kernels::update_triplet_weights const& /* kernel */,
        TVec const& blockThreadExtent, TVec const& /* threadElemExtent */,
        traccc::seedfilter_config config, TArgs const&... /* args */
        ) -> std::size_t {
        return static_cast<std::size_t>(config.compatSeedLimit *
                                        blockThreadExtent.prod()) *
               sizeof(traccc::scalar);
    }
};

template <typename TAcc>
struct BlockSharedMemDynSizeBytes<traccc::alpaka::kernels::select_seeds, TAcc> {
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
        traccc::alpaka::kernels::select_seeds const& /* kernel */,
        TVec const& blockThreadExtent, TVec const& /* threadElemExtent */,
        traccc::seedfinder_config config, TArgs const&... /* args */
        ) -> std::size_t {
        return static_cast<std::size_t>(config.maxSeedsPerSpM *
                                        blockThreadExtent.prod()) *
               sizeof(traccc::device::device_triplet);
    }
};

}  // namespace alpaka::trait
