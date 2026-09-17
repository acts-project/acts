// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/update_tip_length_buffer.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::update_tip_length_buffer
__global__ void update_tip_length_buffer(
    const device::update_tip_length_buffer_payload payload) {
  device::update_tip_length_buffer(details::global_index1(), payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::update_tip_length_buffer_kernel(
    unsigned int n_threads,
    const device::update_tip_length_buffer_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = 512;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  hipLaunchKernelGGL(kernels::update_tip_length_buffer, deviceBlocks,
                     deviceThreads, 0, details::get_stream(stream()), payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
