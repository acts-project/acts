/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/clusterization/clusterization_algorithm.hpp"

#include "../utils/barrier.hpp"
#include "../utils/get_queue.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s)
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/clusterization/device/ccl_kernel.hpp"
#include "traccc/clusterization/device/reify_cluster_data.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/relations.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::alpaka {
namespace kernels {

/// Alpaka kernel for running @c traccc::device::ccl_kernel
struct ccl_kernel {

    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const clustering_config cfg,
        const edm::silicon_cell_collection::const_view cells_view,
        const detector_design_description::const_view det_descr_view,
        const detector_conditions_description::const_view det_cond_view,
        vecmem::data::vector_view<device::details::fallback_index_t>
            f_backup_view,
        vecmem::data::vector_view<device::details::fallback_index_t>
            gf_backup_view,
        vecmem::data::vector_view<unsigned char> adjc_backup_view,
        vecmem::data::vector_view<device::details::fallback_index_t>
            adjv_backup_view,
        uint32_t* backup_mutex_ptr,
        vecmem::data::vector_view<unsigned int> disjoint_set_view,
        vecmem::data::vector_view<unsigned int> cluster_size_view,
        edm::measurement_collection::view measurements_view) const {

        details::thread_id1 thread_id(acc);

        auto& partition_start =
            ::alpaka::declareSharedVar<std::size_t, __COUNTER__>(acc);
        auto& partition_end =
            ::alpaka::declareSharedVar<std::size_t, __COUNTER__>(acc);
        auto& outi = ::alpaka::declareSharedVar<std::size_t, __COUNTER__>(acc);

        device::details::index_t* const shared_v =
            ::alpaka::getDynSharedMem<device::details::index_t>(acc);
        vecmem::data::vector_view<device::details::index_t> f_view{
            cfg.max_partition_size(), shared_v};
        vecmem::data::vector_view<device::details::index_t> gf_view{
            cfg.max_partition_size(), shared_v + cfg.max_partition_size()};

        vecmem::device_atomic_ref<uint32_t> backup_mutex(*backup_mutex_ptr);

        alpaka::barrier<TAcc> barry_r(&acc);

        device::ccl_kernel(cfg, thread_id, cells_view, det_descr_view,
                           det_cond_view, partition_start, partition_end, outi,
                           f_view, gf_view, f_backup_view, gf_backup_view,
                           adjc_backup_view, adjv_backup_view, backup_mutex,
                           disjoint_set_view, cluster_size_view, barry_r,
                           measurements_view);
    }

};  // struct ccl_kernel

/// Alpaka kernel for running @c traccc::device::reify_cluster_data
struct reify_cluster_data {

    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        vecmem::data::vector_view<const unsigned int> disjoint_set_view,
        vecmem::data::vector_view<const unsigned int> permutation_map_view,
        traccc::edm::silicon_cluster_collection::view cluster_view) const {

        device::reify_cluster_data(details::thread_id1{acc}.getGlobalThreadId(),
                                   disjoint_set_view, permutation_map_view,
                                   cluster_view);
    }

};  // struct reify_cluster_data

}  // namespace kernels

clusterization_algorithm::clusterization_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    alpaka::queue& q, const config_type& config,
    std::unique_ptr<const Logger> logger)
    : device::clusterization_algorithm(mr, copy, config, std::move(logger)),
      alpaka::algorithm_base(q) {}

bool clusterization_algorithm::input_is_contiguous(
    const edm::silicon_cell_collection::const_view&) const {

    // TODO: implement sanity checks for the input data in Alpaka
    return true;
}

bool clusterization_algorithm::input_is_sorted(
    const edm::silicon_cell_collection::const_view&) const {

    // TODO: implement sanity checks for the input data in Alpaka
    return true;
}

void clusterization_algorithm::sort_cells_kernel(
    const unsigned int, const edm::silicon_cell_collection::const_view&,
    edm::silicon_cell_collection::view&,
    vecmem::data::vector_view<unsigned int>&, const bool) const {

    throw std::runtime_error(
        "Cell sorting is not yet implemented for the Alpaka backend.");
}

void clusterization_algorithm::ccl_kernel(
    const ccl_kernel_payload& payload) const {

    Idx num_blocks =
        (payload.n_cells + (payload.config.target_partition_size()) - 1) /
        payload.config.target_partition_size();
    static_assert(::alpaka::isMultiThreadAcc<Acc>,
                  "Clustering algorithm must be compiled for an accelerator "
                  "with support for multi-thread blocks.");
    auto workDiv =
        makeWorkDiv<Acc>(num_blocks, payload.config.threads_per_partition);
    ::alpaka::exec<Acc>(
        details::get_queue(queue()), workDiv, kernels::ccl_kernel{},
        payload.config, payload.cells, payload.det_descr, payload.det_cond,
        payload.f_backup, payload.gf_backup, payload.adjc_backup,
        payload.adjv_backup, payload.backup_mutex, payload.disjoint_set,
        payload.cluster_sizes, payload.measurements);
}

void clusterization_algorithm::cluster_maker_kernel(
    unsigned int num_cells,
    const vecmem::data::vector_view<unsigned int>& disjoint_set,
    edm::silicon_cluster_collection::view& cluster_data,
    const vecmem::data::vector_view<const unsigned int>& permutation_map_view)
    const {

    const unsigned int num_threads = warp_size() * 16u;
    const unsigned int num_blocks = (num_cells + num_threads - 1) / num_threads;
    static_assert(::alpaka::isMultiThreadAcc<Acc>,
                  "Clustering algorithm must be compiled for an accelerator "
                  "with support for multi-thread blocks.");
    auto workDiv = makeWorkDiv<Acc>(num_blocks, num_threads);
    ::alpaka::exec<Acc>(details::get_queue(queue()), workDiv,
                        kernels::reify_cluster_data{}, disjoint_set,
                        permutation_map_view, cluster_data);
    // The base class destroys the input buffers right after this call, so
    // the kernel must finish before returning.
    queue().synchronize();
}

}  // namespace traccc::alpaka

// Define the required trait needed for Dynamic shared memory allocation.
namespace alpaka::trait {

template <typename TAcc>
struct BlockSharedMemDynSizeBytes<traccc::alpaka::kernels::ccl_kernel, TAcc> {
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
        traccc::alpaka::kernels::ccl_kernel const& /* kernel */,
        TVec const& /* blockThreadExtent */, TVec const& /* threadElemExtent */,
        const traccc::clustering_config config, TArgs const&... /* args */
        ) -> std::size_t {
        return static_cast<std::size_t>(
            2 * config.max_partition_size() *
            sizeof(traccc::device::details::index_t));
    }
};

}  // namespace alpaka::trait
