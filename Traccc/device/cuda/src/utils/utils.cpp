/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "utils.hpp"

#include "cuda_error_handling.hpp"

namespace traccc::cuda::details {

unsigned int get_warp_size(int device) {

    int warp_size = 0;
    TRACCC_CUDA_ERROR_CHECK(
        cudaDeviceGetAttribute(&warp_size, cudaDevAttrWarpSize, device));
    return static_cast<unsigned int>(warp_size);
}

cudaStream_t get_stream(const stream& stream) {

    return reinterpret_cast<cudaStream_t>(stream.cudaStream());
}

}  // namespace traccc::cuda::details
