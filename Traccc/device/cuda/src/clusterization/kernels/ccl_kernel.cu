/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// CUDA Library include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/cuda_error_handling.hpp"
#include "../../utils/thread_id.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"

// Project include(s)
#include "traccc/clusterization/device/ccl_kernel.hpp"

namespace traccc::cuda::kernels {

/// CUDA kernel for running @c traccc::device::ccl_kernel
__global__ void ccl_kernel(
    const clustering_config cfg,
    const edm::silicon_cell_collection::const_view cells_view,
    const detector_design_description::const_view det_desc_view,
    const detector_conditions_description::const_view det_cond_view,
    edm::measurement_collection::view measurements_view,
    vecmem::data::vector_view<device::details::fallback_index_t> f_backup_view,
    vecmem::data::vector_view<device::details::fallback_index_t> gf_backup_view,
    vecmem::data::vector_view<unsigned char> adjc_backup_view,
    vecmem::data::vector_view<device::details::fallback_index_t>
        adjv_backup_view,
    unsigned int* backup_mutex_ptr,
    vecmem::data::vector_view<unsigned int> disjoint_set_view,
    vecmem::data::vector_view<unsigned int> cluster_size_view) {
  __shared__ std::size_t partition_start, partition_end;
  __shared__ std::size_t outi;
  extern __shared__ device::details::index_t shared_v[];
  vecmem::device_atomic_ref<unsigned int> backup_mutex(*backup_mutex_ptr);

  using vector_size_t =
      vecmem::data::vector_view<device::details::index_t>::size_type;

  vecmem::data::vector_view<device::details::index_t> f_view{
      static_cast<vector_size_t>(cfg.max_partition_size()), shared_v};
  vecmem::data::vector_view<device::details::index_t> gf_view{
      static_cast<vector_size_t>(cfg.max_partition_size()),
      shared_v + cfg.max_partition_size()};
  traccc::cuda::barrier barry_r;
  const details::thread_id1 thread_id;

  device::ccl_kernel(cfg, thread_id, cells_view, det_desc_view, det_cond_view,
                     partition_start, partition_end, outi, f_view, gf_view,
                     f_backup_view, gf_backup_view, adjc_backup_view,
                     adjv_backup_view, backup_mutex, disjoint_set_view,
                     cluster_size_view, barry_r, measurements_view);
}
}  // namespace traccc::cuda::kernels
