// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "../create_device_detector.cuh"

// System include(s).
#include <new>

namespace traccc::cuda {
namespace kernels {

template <detray::concepts::detector detector_t>
__global__ void create_device_detector(
    const detray::detector_view_t<detector_t> in, detector_t* out) {
  unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (thread_id == 0) {
    new (out) detector_t(in);
  }
}

}  // namespace kernels

template <detray::concepts::detector detector_t>
void create_device_detector(const cudaStream_t& stream,
                            const detray::detector_view_t<detector_t> in,
                            detector_t* out) {
  kernels::create_device_detector<detector_t><<<1, 1, 0, stream>>>(in, out);
}
}  // namespace traccc::cuda
