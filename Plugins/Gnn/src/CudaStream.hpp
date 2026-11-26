// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/Tensor.hpp"

#include <cassert>

#include <cuda_runtime_api.h>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

namespace ActsPlugins::detail {

/// Get the `cudaStream_t` of an execution context, which must carry a stream
inline cudaStream_t getCudaStream(const ExecutionContext &ctx) {
  assert(ctx.stream != nullptr);
  return static_cast<cudaStream_t>(ctx.stream->stream());
}

}  // namespace ActsPlugins::detail
