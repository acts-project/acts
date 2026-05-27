/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// CUDA Library include(s).
#include "../sanity/contiguous_on.cuh"
#include "../sanity/ordered_on.cuh"
#include "../utils/cuda_error_handling.hpp"
#include "./kernels/ccl_kernel.cuh"
#include "./kernels/reify_cluster_data.cuh"
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/relations.hpp"

// Vecmem include(s).
#include <cstring>
#include <vecmem/utils/copy.hpp>

namespace traccc::cuda {

clusterization_algorithm::clusterization_algorithm(
    const traccc::memory_resource& mr, vecmem::copy& copy, cuda::stream& str,
    const config_type& config, std::unique_ptr<const Logger> logger)
    : device::clusterization_algorithm(mr, copy, config, std::move(logger)),
      cuda::algorithm_base(str) {}

bool clusterization_algorithm::input_is_valid(
    const edm::silicon_cell_collection::const_view& cells) const {

    return (is_contiguous_on<edm::silicon_cell_collection::const_device>(
                cell_module_projection(), mr().main, copy(), stream(), cells) &&
            is_ordered_on<edm::silicon_cell_collection::const_device>(
                channel0_major_cell_order_relation(), mr().main, copy(),
                stream(), cells));
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
}

void clusterization_algorithm::cluster_maker_kernel(
    unsigned int num_cells,
    const vecmem::data::vector_view<unsigned int>& disjoint_set,
    edm::silicon_cluster_collection::view& cluster_data) const {

    const unsigned int num_threads = warp_size() * 16u;
    const unsigned int num_blocks = (num_cells + num_threads - 1) / num_threads;
    kernels::reify_cluster_data<<<num_blocks, num_threads, 0,
                                  details::get_stream(stream())>>>(
        disjoint_set, cluster_data);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
