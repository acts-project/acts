/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "get_device.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda::details {

int get_device() {

    int d = -1;
    cudaGetDevice(&d);
    return d;
}

}  // namespace traccc::cuda::details
