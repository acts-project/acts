/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cuda_runtime.h>

#include "traccc/edm/track_container.hpp"

namespace traccc::cuda {
void fit_prelude(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream,
    vecmem::data::vector_view<const unsigned int> param_ids_view,
    edm::track_container<default_algebra>::const_view track_candidates_view,
    edm::track_container<default_algebra>::view tracks_view,
    vecmem::data::vector_view<unsigned int> param_liveness_view);
}
