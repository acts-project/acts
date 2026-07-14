/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// CUDA Library include(s).
#include "../sanity/contiguous_on.cuh"
#include "../sanity/ordered_on.cuh"
#include "../utils/barrier.hpp"
#include "../utils/cuda_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "./kernels/ccl_kernel.cuh"
#include "./kernels/reify_cluster_data.cuh"
#include "./kernels/sort_cells.cuh"
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/relations.hpp"

// Vecmem include(s).
#include <cstring>
#include <limits>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/utils/copy.hpp>

namespace traccc::cuda {

clusterization_algorithm::clusterization_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    const stream_wrapper& str, const config_type& config,
    std::unique_ptr<const Logger> logger)
    : device::clusterization_algorithm(mr, copy, config, std::move(logger)),
      cuda::algorithm_base(str) {}

bool clusterization_algorithm::input_is_contiguous(
    const edm::silicon_cell_collection::const_view& cells) const {

    return is_contiguous_on<edm::silicon_cell_collection::const_device>(
        cell_module_projection(), mr().main, copy(), stream(), cells);
}

bool clusterization_algorithm::input_is_sorted(
    const edm::silicon_cell_collection::const_view& cells) const {

    return is_ordered_on<edm::silicon_cell_collection::const_device>(
        channel0_major_cell_order_relation(), mr().main, copy(), stream(),
        cells);
}

void clusterization_algorithm::sort_cells_kernel(
    const unsigned int num_cells,
    const edm::silicon_cell_collection::const_view& cells,
    edm::silicon_cell_collection::view& new_cells,
    vecmem::data::vector_view<unsigned int>& permutation_map_view,
    const bool write_permutation) const {

    static constexpr std::size_t SORT_THREADS_PER_BLOCK = 512;

    const unsigned blockSize = SORT_THREADS_PER_BLOCK;
    const unsigned int cellsPerThread = 4;
    const unsigned int numBlocks =
        (num_cells + (blockSize * cellsPerThread) - 1) /
        (blockSize * cellsPerThread);

    kernels::sort_cells<SORT_THREADS_PER_BLOCK, 8>
        <<<numBlocks, blockSize, 0, details::get_stream(stream())>>>(
            cellsPerThread, num_cells, cells, new_cells, permutation_map_view,
            write_permutation);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void clusterization_algorithm::ccl_kernel(
    const ccl_kernel_payload& payload) const {

    const unsigned int num_blocks =
        (payload.n_cells + (payload.config.target_partition_size()) - 1) /
        payload.config.target_partition_size();
    kernels::ccl_kernel<<<num_blocks, payload.config.threads_per_partition,
                          2 * payload.config.max_partition_size() *
                              sizeof(device::details::index_t),
                          details::get_stream(stream())>>>(
        payload.config, payload.cells, payload.det_descr, payload.det_cond,
        payload.measurements, payload.f_backup, payload.gf_backup,
        payload.adjc_backup, payload.adjv_backup, payload.backup_mutex,
        payload.disjoint_set, payload.cluster_sizes);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    // With sorting enabled, the kernel reads from sorted cell and permutation
    // map buffers that the base class destroys without any further
    // synchronization, so the kernel must finish before returning.
    if (payload.config.sort_cells) {
        stream().synchronize();
    }
}

void clusterization_algorithm::cluster_maker_kernel(
    unsigned int num_cells,
    const vecmem::data::vector_view<unsigned int>& disjoint_set,
    edm::silicon_cluster_collection::view& cluster_data,
    const vecmem::data::vector_view<const unsigned int>& permutation_map_view)
    const {

    const unsigned int num_threads = warp_size() * 16u;
    const unsigned int num_blocks = (num_cells + num_threads - 1) / num_threads;
    kernels::reify_cluster_data<<<num_blocks, num_threads, 0,
                                  details::get_stream(stream())>>>(
        disjoint_set, permutation_map_view, cluster_data);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    // The base class destroys the input buffers right after this call, so
    // the kernel must finish before returning.
    stream().synchronize();
}

}  // namespace traccc::cuda
