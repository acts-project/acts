/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/device/count_doublets.hpp"
#include "traccc/seeding/device/count_grid_capacities.hpp"
#include "traccc/seeding/device/count_triplets.hpp"
#include "traccc/seeding/device/find_doublets.hpp"
#include "traccc/seeding/device/find_triplets.hpp"
#include "traccc/seeding/device/populate_grid.hpp"
#include "traccc/seeding/device/reduce_triplet_counts.hpp"
#include "traccc/seeding/device/select_seeds.hpp"
#include "traccc/seeding/device/update_triplet_weights.hpp"

namespace traccc::cuda {
namespace kernels {

/// CUDA kernel for running @c traccc::device::count_grid_capacities
__global__ void count_grid_capacities(
    seedfinder_config config,
    traccc::details::spacepoint_grid_types::host::axis_p0_type phi_axis,
    traccc::details::spacepoint_grid_types::host::axis_p1_type z_axis,
    edm::spacepoint_collection::const_view spacepoints,
    vecmem::data::vector_view<unsigned int> grid_capacities) {

    device::count_grid_capacities(details::global_index1(), config, phi_axis,
                                  z_axis, spacepoints, grid_capacities);
}

/// CUDA kernel for running @c traccc::device::populate_grid
__global__ void populate_grid(
    seedfinder_config config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::view grid,
    vecmem::data::vector_view<device::prefix_sum_element_t> grid_prefix_sum) {

    device::populate_grid(details::global_index1(), config, spacepoints, grid,
                          grid_prefix_sum);
}

/// CUDA kernel for running @c traccc::device::count_doublets
__global__ void count_doublets(
    seedfinder_config config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::const_view sp_grid,
    vecmem::data::vector_view<const device::prefix_sum_element_t> sp_prefix_sum,
    device::doublet_counter_collection_types::view doublet_counter,
    unsigned int& nMidBot, unsigned int& nMidTop) {

    device::count_doublets(details::global_index1(), config, spacepoints,
                           sp_grid, sp_prefix_sum, doublet_counter, nMidBot,
                           nMidTop);
}

/// CUDA kernel for running @c traccc::device::find_doublets
__global__ void find_doublets(
    seedfinder_config config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::const_view sp_grid,
    device::doublet_counter_collection_types::const_view doublet_counter,
    device::device_doublet_collection_types::view mb_doublets,
    device::device_doublet_collection_types::view mt_doublets) {

    device::find_doublets(details::global_index1(), config, spacepoints,
                          sp_grid, doublet_counter, mb_doublets, mt_doublets);
}

/// CUDA kernel for running @c traccc::device::count_triplets
__global__ void count_triplets(
    seedfinder_config config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::const_view sp_grid,
    device::doublet_counter_collection_types::const_view doublet_counter,
    device::device_doublet_collection_types::const_view mb_doublets,
    device::device_doublet_collection_types::const_view mt_doublets,
    device::triplet_counter_spM_collection_types::view spM_counter,
    device::triplet_counter_collection_types::view midBot_counter) {

    device::count_triplets(details::global_index1(), config, spacepoints,
                           sp_grid, doublet_counter, mb_doublets, mt_doublets,
                           spM_counter, midBot_counter);
}

/// CUDA kernel for running @c traccc::device::reduce_triplet_counts
__global__ void reduce_triplet_counts(
    device::doublet_counter_collection_types::const_view doublet_counter,
    device::triplet_counter_spM_collection_types::view spM_counter,
    unsigned int& num_triplets) {

    device::reduce_triplet_counts(details::global_index1(), doublet_counter,
                                  spM_counter, num_triplets);
}

/// CUDA kernel for running @c traccc::device::find_triplets
__global__ void find_triplets(
    seedfinder_config config, seedfilter_config filter_config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::const_view sp_grid,
    device::doublet_counter_collection_types::const_view doublet_counter,
    device::device_doublet_collection_types::const_view mt_doublets,
    device::triplet_counter_spM_collection_types::const_view spM_tc,
    device::triplet_counter_collection_types::const_view midBot_tc,
    device::device_triplet_collection_types::view triplet_view) {

    device::find_triplets(details::global_index1(), config, filter_config,
                          spacepoints, sp_grid, doublet_counter, mt_doublets,
                          spM_tc, midBot_tc, triplet_view);
}

/// CUDA kernel for running @c traccc::device::update_triplet_weights
__global__ void update_triplet_weights(
    seedfilter_config filter_config,
    edm::spacepoint_collection::const_view spacepoints,
    device::triplet_counter_spM_collection_types::const_view spM_tc,
    device::triplet_counter_collection_types::const_view midBot_tc,
    device::device_triplet_collection_types::view triplet_view) {

    // Array for temporary storage of quality parameters for comparing triplets
    // within weight updating kernel
    extern __shared__ scalar data[];
    // Each thread uses compatSeedLimit elements of the array
    scalar* dataPos = &data[threadIdx.x * filter_config.compatSeedLimit];

    device::update_triplet_weights(details::global_index1(), filter_config,
                                   spacepoints, spM_tc, midBot_tc, dataPos,
                                   triplet_view);
}

/// CUDA kernel for running @c traccc::device::select_seeds
__global__ void select_seeds(
    seedfinder_config finder_config, seedfilter_config filter_config,
    edm::spacepoint_collection::const_view spacepoints,
    traccc::details::spacepoint_grid_types::const_view sp_view,
    device::triplet_counter_spM_collection_types::const_view spM_tc,
    device::triplet_counter_collection_types::const_view midBot_tc,
    device::device_triplet_collection_types::const_view triplet_view,
    edm::seed_collection::view seed_view) {

    // Array for temporary storage of triplets for comparing within seed
    // selecting kernel
    extern __shared__ device::device_triplet data2[];
    // Each thread uses max_triplets_per_spM elements of the array
    device::device_triplet* dataPos =
        &data2[threadIdx.x * finder_config.maxSeedsPerSpM];

    device::select_seeds(details::global_index1(), finder_config, filter_config,
                         spacepoints, sp_view, spM_tc, midBot_tc, triplet_view,
                         dataPos, seed_view);
}

}  // namespace kernels

triplet_seeding_algorithm::triplet_seeding_algorithm(
    const seedfinder_config& finder_config,
    const spacepoint_grid_config& grid_config,
    const seedfilter_config& filter_config, const traccc::memory_resource& mr,
    vecmem::copy& copy, cuda::stream& str, std::unique_ptr<const Logger> logger)
    : device::triplet_seeding_algorithm(finder_config, grid_config,
                                        filter_config, mr, copy,
                                        std::move(logger)),
      cuda::algorithm_base{str} {}

void triplet_seeding_algorithm::count_grid_capacities_kernel(
    const count_grid_capacities_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    kernels::count_grid_capacities<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(
        payload.config, payload.phi_axis, payload.z_axis, payload.spacepoints,
        payload.grid_capacities);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::populate_grid_kernel(
    const populate_grid_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    kernels::populate_grid<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(
        payload.config, payload.spacepoints, payload.grid,
        payload.grid_prefix_sum);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::count_doublets_kernel(
    const count_doublets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_spacepoints + n_threads - 1) / n_threads;
    kernels::count_doublets<<<n_blocks, n_threads, 0,
                              details::get_stream(stream())>>>(
        payload.config, payload.spacepoints, payload.grid,
        payload.grid_prefix_sum, payload.doublet_counter, payload.nMidBot,
        payload.nMidTop);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::find_doublets_kernel(
    const find_doublets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    kernels::find_doublets<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(
        payload.config, payload.spacepoints, payload.grid,
        payload.doublet_counter, payload.mb_doublets, payload.mt_doublets);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::count_triplets_kernel(
    const count_triplets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks = (payload.nMidBot + n_threads - 1) / n_threads;
    kernels::count_triplets<<<n_blocks, n_threads, 0,
                              details::get_stream(stream())>>>(
        payload.config, payload.spacepoints, payload.grid,
        payload.doublet_counter, payload.mb_doublets, payload.mt_doublets,
        payload.spM_counter, payload.midBot_counter);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::triplet_counts_reduction_kernel(
    const triplet_counts_reduction_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    kernels::reduce_triplet_counts<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(
        payload.doublet_counter, payload.spM_counter, payload.nTriplets);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::find_triplets_kernel(
    const find_triplets_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks = (payload.nMidBot + n_threads - 1) / n_threads;
    kernels::find_triplets<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(
        payload.finding_config, payload.filter_config, payload.spacepoints,
        payload.grid, payload.doublet_counter, payload.mt_doublets,
        payload.spM_tc, payload.midBot_tc, payload.triplets);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::update_triplet_weights_kernel(
    const update_triplet_weights_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_triplets + n_threads - 1) / n_threads;
    kernels::update_triplet_weights<<<
        n_blocks, n_threads,
        sizeof(scalar) * payload.config.compatSeedLimit * n_threads,
        details::get_stream(stream())>>>(payload.config, payload.spacepoints,
                                         payload.spM_tc, payload.midBot_tc,
                                         payload.triplets);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void triplet_seeding_algorithm::select_seeds_kernel(
    const select_seeds_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 2;
    const unsigned int n_blocks =
        (payload.n_doublets + n_threads - 1) / n_threads;
    kernels::
        select_seeds<<<n_blocks, n_threads,
                       sizeof(device::device_triplet) *
                           payload.finder_config.maxSeedsPerSpM * n_threads,
                       details::get_stream(stream())>>>(
            payload.finder_config, payload.filter_config, payload.spacepoints,
            payload.grid, payload.spM_tc, payload.midBot_tc, payload.triplets,
            payload.seeds);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
