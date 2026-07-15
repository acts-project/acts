/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "./kernels/ccl_kernel.hpp"
#include "./kernels/reify_cluster_data.hpp"
#include "traccc/hip/clusterization/clusterization_algorithm.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <stdexcept>

namespace traccc::hip {

clusterization_algorithm::clusterization_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    hip::stream& str, const config_type& config,
    std::unique_ptr<const Logger> logger)
    : device::clusterization_algorithm(mr, copy, config, std::move(logger)),
      hip::algorithm_base(str) {}

bool clusterization_algorithm::input_is_contiguous(
    const edm::silicon_cell_collection::const_view& /*cells*/) const {

    /// TODO: implement input checks
    return true;
}

bool clusterization_algorithm::input_is_sorted(
    const edm::silicon_cell_collection::const_view& /*cells*/) const {

    /// TODO: implement input checks
    return true;
}

void clusterization_algorithm::sort_cells_kernel(
    const unsigned int, const edm::silicon_cell_collection::const_view&,
    edm::silicon_cell_collection::view&,
    vecmem::data::vector_view<unsigned int>&, const bool) const {

    throw std::runtime_error(
        "Cell sorting is not yet implemented for the HIP backend.");
}

void clusterization_algorithm::ccl_kernel(
    const ccl_kernel_payload& payload) const {

    const unsigned int num_blocks =
        (payload.n_cells + (payload.config.target_partition_size()) - 1) /
        payload.config.target_partition_size();
    hipLaunchKernelGGL(kernels::ccl_kernel, dim3{num_blocks},
                       dim3{payload.config.threads_per_partition},
                       2 * payload.config.max_partition_size() *
                           sizeof(device::details::index_t),
                       details::get_stream(stream()), payload.config,
                       payload.cells, payload.det_descr, payload.det_cond,
                       payload.measurements, payload.f_backup,
                       payload.gf_backup, payload.adjc_backup,
                       payload.adjv_backup, payload.backup_mutex,
                       payload.disjoint_set, payload.cluster_sizes);
    TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

void clusterization_algorithm::cluster_maker_kernel(
    unsigned int num_cells,
    const vecmem::data::vector_view<unsigned int>& disjoint_set,
    edm::silicon_cluster_collection::view& cluster_data,
    const vecmem::data::vector_view<const unsigned int>& permutation_map_view)
    const {

    const unsigned int num_threads = warp_size() * 16u;
    const unsigned int num_blocks = (num_cells + num_threads - 1) / num_threads;
    hipLaunchKernelGGL(kernels::reify_cluster_data, num_blocks, num_threads, 0,
                       details::get_stream(stream()), disjoint_set,
                       permutation_map_view, cluster_data);
    TRACCC_HIP_ERROR_CHECK(hipGetLastError());
    // The base class destroys the input buffers right after this call, so
    // the kernel must finish before returning.
    stream().synchronize();
}

}  // namespace traccc::hip
