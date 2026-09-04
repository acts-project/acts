/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::cuda::kernels {

/// CUDA kernel for running @c traccc::device::ccl_kernel
__global__ void ccl_kernel(
    const clustering_config cfg,
    const edm::silicon_cell_collection::const_view cells_view,
    const detector_design_description::const_view det_descr_view,
    const detector_conditions_description::const_view det_cond_view,
    edm::measurement_collection::view measurements_view,
    vecmem::data::vector_view<device::details::fallback_index_t> f_backup_view,
    vecmem::data::vector_view<device::details::fallback_index_t> gf_backup_view,
    vecmem::data::vector_view<unsigned char> adjc_backup_view,
    vecmem::data::vector_view<device::details::fallback_index_t>
        adjv_backup_view,
    unsigned int* backup_mutex_ptr,
    vecmem::data::vector_view<unsigned int> disjoint_set_view,
    vecmem::data::vector_view<unsigned int> cluster_size_view);

}  // namespace traccc::cuda::kernels
