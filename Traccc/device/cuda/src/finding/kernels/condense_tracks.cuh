/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/condense_tracks.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

void condense_tracks(const dim3& grid_size, const dim3& block_size,
                     std::size_t shared_mem_size, const cudaStream_t& stream,
                     const device::condense_tracks_payload& payload);

}  // namespace traccc::cuda
