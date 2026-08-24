/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../create_device_detector.cuh"

// System include(s).
#include <new>

namespace traccc::cuda {
namespace kernels {

template <typename detector_t>
__global__ void create_device_detector(typename detector_t::const_view_type in,
                                       detector_t* out) {
  unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (thread_id == 0) {
    new (out) detector_t(in);
  }
}

}  // namespace kernels

template <typename detector_t>
void create_device_detector(const cudaStream_t& stream,
                            typename detector_t::const_view_type in,
                            detector_t* out) {
  kernels::create_device_detector<detector_t><<<1, 1, 0, stream>>>(in, out);
}
}  // namespace traccc::cuda
