// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts::detail {

template <typename T>
__global__ void iota(std::size_t size, T *array) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) {
    return;
  }
  array[i] = i;
}

}  // namespace Acts::detail
