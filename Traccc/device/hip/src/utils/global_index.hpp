// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/device/global_index.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip::details {

/// Function creating a global index in a 1D HIP kernel
__device__ inline device::global_index_t global_index1() {
  return blockIdx.x * blockDim.x + threadIdx.x;
}

}  // namespace traccc::hip::details
