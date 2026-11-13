// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/global_index.hpp"
#include "permute_tracks.cuh"

// Project include(s).
#include "traccc/finding/device/permute_tracks.hpp"

namespace traccc::cuda {
namespace kernels {

__global__ void permute_tracks(
    const __grid_constant__ device::permute_tracks_payload payload) {
  device::permute_tracks(details::global_index1(), payload);
}

}  // namespace kernels

void permute_tracks(const dim3& grid_size, const dim3& block_size,
                    std::size_t shared_mem_size, const cudaStream_t& stream,
                    const device::permute_tracks_payload& payload) {
  kernels::permute_tracks<<<grid_size, block_size, shared_mem_size, stream>>>(
      payload);
}
}  // namespace traccc::cuda
