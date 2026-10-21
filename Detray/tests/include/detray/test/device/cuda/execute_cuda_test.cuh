// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/cuda_definitions.hpp"

// System include(s).
#include <cstddef>

namespace detray::test::cuda {

template <class functor_t, typename... Args>
__global__ void cuda_test_kernel(std::size_t array_sizes, Args... args) {
  // Find the current index that we need to process.
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= array_sizes) {
    return;
  }

  // Execute the test functor for this index.
  functor_t()(i, std::forward<Args>(args)...);
}

/// Execute a test functor on a device, on @c array_sizes threads
template <class functor_t, class... Args>
void execute_cuda_test(std::size_t array_sizes, Args... args) {
  // Number of threads per execution block. Less than 1024 to make debug tests
  // possible.
  const int n_threads_per_block{std::min(256, static_cast<int>(array_sizes))};
  const int n_blocks{(static_cast<int>(array_sizes) + n_threads_per_block - 1) /
                     n_threads_per_block};

  // Launch the test on the device.
  cuda_test_kernel<functor_t><<<n_blocks, n_threads_per_block>>>(
      array_sizes, std::forward<Args>(args)...);

  // Check whether it succeeded to run.
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray::test::cuda
