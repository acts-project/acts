/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/device/global_index.hpp"

namespace traccc::hip::details {

/// Function creating a global index in a 1D HIP kernel
__device__ inline device::global_index_t global_index1() {
  return blockIdx.x * blockDim.x + threadIdx.x;
}

}  // namespace traccc::hip::details
