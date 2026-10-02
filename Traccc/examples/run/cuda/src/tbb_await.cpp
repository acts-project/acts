// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/cuda/tbb_await.hpp"

#include "cuda_error_check.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

// TBB include(s).
#include <tbb/task.h>

namespace traccc::cuda {

namespace {
void CUDART_CB suspend_stream_callback(void* tag) {
  tbb::task::resume(*static_cast<tbb::task::suspend_point*>(tag));
}
}  // namespace

void tbb_await_callback(vecmem::abstract_event& event,
                        const stream_wrapper& stream) {
  event.ignore();  // ignore the event, as it is not needed for resumption
  cudaError_t err = cudaSuccess;
  auto suspend_point =
      tbb::task::suspend_point{};  // suspension point address must remain valid
                                   // when resumption callback is called
  tbb::task::suspend([&err, &stream, &suspend_point](auto tag) {
    suspend_point = tag;
    auto cuda_stream = reinterpret_cast<cudaStream_t>(stream.cudaStream());
    err = cudaLaunchHostFunc(cuda_stream, suspend_stream_callback,
                             &suspend_point);
    // resume immediately if the callback could not be registered
    if (err != cudaSuccess) {
      tbb::task::resume(suspend_point);
    }
  });
  CUDA_ERROR_CHECK(err);
}
}  // namespace traccc::cuda
