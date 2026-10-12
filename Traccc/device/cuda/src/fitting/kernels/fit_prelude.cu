// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/cuda_error_handling.hpp"
#include "../../utils/global_index.hpp"
#include "./fit_prelude.hpp"

// Project include(s).
#include "traccc/fitting/device/fit_prelude.hpp"

namespace traccc::cuda {
namespace kernels {

__global__ void fit_prelude(const device::fit_prelude_payload payload) {
  device::fit_prelude(details::global_index1(), payload);
}

}  // namespace kernels

void fit_prelude(const dim3& grid_size, const dim3& block_size,
                 std::size_t shared_mem_size, const cudaStream_t& stream,
                 const device::fit_prelude_payload& payload) {
  kernels::fit_prelude<<<grid_size, block_size, shared_mem_size, stream>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
