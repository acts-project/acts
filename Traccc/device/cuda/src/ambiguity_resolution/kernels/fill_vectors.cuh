/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"
#include "traccc/ambiguity_resolution/device/fill_vectors.hpp"

namespace traccc::cuda::kernels {

__global__ void fill_vectors(const ambiguity_resolution_config cfg,
                             device::fill_vectors_payload payload);
}
