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
#include "traccc/finding/device/remove_duplicates.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::remove_duplicates
__global__ void remove_duplicates(
    const finding_config cfg, const device::remove_duplicates_payload payload) {
  device::remove_duplicates(details::global_index1(), cfg, payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::remove_duplicates_kernel(
    unsigned int n_threads, const finding_config& config,
    const device::remove_duplicates_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 8;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  hipLaunchKernelGGL(kernels::remove_duplicates, deviceBlocks, deviceThreads, 0,
                     details::get_stream(stream()), config, payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
