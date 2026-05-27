/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "../../utils/global_index.hpp"
#include "./fit_prelude.hpp"
#include "traccc/fitting/device/fit_prelude.hpp"

namespace traccc::cuda {
namespace kernels {
__global__ void fit_prelude(
    vecmem::data::vector_view<const unsigned int> param_ids_view,
    edm::track_container<default_algebra>::const_view track_candidates_view,
    edm::track_container<default_algebra>::view tracks_view,
    vecmem::data::vector_view<unsigned int> param_liveness_view) {
    device::fit_prelude<default_algebra>(details::global_index1(),
                                         param_ids_view, track_candidates_view,
                                         tracks_view, param_liveness_view);
}
}  // namespace kernels

void fit_prelude(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream,
    vecmem::data::vector_view<const unsigned int> param_ids_view,
    edm::track_container<default_algebra>::const_view track_candidates_view,
    edm::track_container<default_algebra>::view tracks_view,
    vecmem::data::vector_view<unsigned int> param_liveness_view) {
    kernels::fit_prelude<<<grid_size, block_size, shared_mem_size, stream>>>(
        param_ids_view, track_candidates_view, tracks_view,
        param_liveness_view);
}
}  // namespace traccc::cuda
