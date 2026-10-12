// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/cuda/utils/stream_wrapper.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda::details {

/// Get the warp size for a given device.
///
/// @param device The device to query.
///
/// @return The warp size for the device.
///
unsigned int get_warp_size(int device);

/// Get concrete @c cudaStream_t object out of our wrapper
cudaStream_t get_stream(const stream_wrapper& str);

}  // namespace traccc::cuda::details
