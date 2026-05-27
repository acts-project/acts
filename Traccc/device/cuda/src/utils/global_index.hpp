/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/device/global_index.hpp"

namespace traccc::cuda::details {

/// Function creating a global index in a 1D CUDA kernel
__device__ inline device::global_index_t global_index1() {

    return blockIdx.x * blockDim.x + threadIdx.x;
}

}  // namespace traccc::cuda::details
