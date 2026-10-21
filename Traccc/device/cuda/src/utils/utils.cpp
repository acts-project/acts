// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

cudaStream_t get_stream(const stream_wrapper& stream) {
  return static_cast<cudaStream_t>(stream.cudaStream());
}

}  // namespace traccc::cuda::details
