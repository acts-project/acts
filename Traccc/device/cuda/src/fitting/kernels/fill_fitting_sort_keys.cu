/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/cuda_error_handling.hpp"
#include "../../utils/global_index.hpp"
#include "fill_fitting_sort_keys.hpp"

// Project include(s).
#include "traccc/fitting/device/fill_fitting_sort_keys.hpp"

namespace traccc::cuda {
namespace kernels {

__global__ void fill_fitting_sort_keys(
    edm::track_collection<default_algebra>::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {

    device::fill_fitting_sort_keys(details::global_index1(),
                                   track_candidates_view, keys_view, ids_view);
}

}  // namespace kernels

void fill_fitting_sort_keys(
    const dim3& grid_size, const dim3& block_size, cudaStream_t stream,
    edm::track_collection<default_algebra>::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {

    kernels::fill_fitting_sort_keys<<<grid_size, block_size, 0, stream>>>(
        track_candidates_view, keys_view, ids_view);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
