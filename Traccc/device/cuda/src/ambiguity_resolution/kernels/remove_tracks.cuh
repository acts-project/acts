/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/device/remove_tracks.hpp"

namespace traccc::cuda::kernels {

__global__ void remove_tracks(device::remove_tracks_payload payload);
}
