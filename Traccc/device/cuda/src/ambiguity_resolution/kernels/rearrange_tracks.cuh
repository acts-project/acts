/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/device/rearrange_tracks.hpp"

namespace traccc::cuda::kernels {

constexpr const int nThreads_per_track = 4;

__global__ void rearrange_tracks(device::rearrange_tracks_payload payload);
}  // namespace traccc::cuda::kernels
