/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/stream.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda::details {

/// Get the warp size for a given device.
///
/// @param device The device to query.
///
/// @return The warp size for the device.
///
unsigned int get_warp_size(int device);

/// Get concrete @c cudaStream_t object out of our wrapper
cudaStream_t get_stream(const stream& str);

}  // namespace traccc::cuda::details
