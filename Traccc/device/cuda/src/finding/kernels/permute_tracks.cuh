// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/finding/device/permute_tracks.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

void permute_tracks(const dim3& grid_size, const dim3& block_size,
                    std::size_t shared_mem_size, const cudaStream_t& stream,
                    const device::permute_tracks_payload& payload);

}  // namespace traccc::cuda
