// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/hip/tbb_await.hpp"

#include "hip_error_check.hpp"

// HIP include(s).
#include <hip/hip_runtime_api.h>

// TBB include(s).
#include <tbb/task.h>

namespace traccc::hip {

namespace {
void suspend_stream_callback(void* tag) {
  tbb::task::resume(*static_cast<tbb::task::suspend_point*>(tag));
}
}  // namespace

void tbb_await_callback(vecmem::abstract_event& event,
                        const stream_wrapper& stream) {
  event.ignore();  // ignore the event, as it is not needed for resumption
  hipError_t err = hipSuccess;
  auto suspend_point =
      tbb::task::suspend_point{};  // suspension point address must remain valid
                                   // when resumption callback is called
  tbb::task::suspend([&err, &stream, &suspend_point](auto tag) {
    suspend_point = tag;
    auto hip_stream = reinterpret_cast<hipStream_t>(stream.hipStream());
    err =
        hipLaunchHostFunc(hip_stream, suspend_stream_callback, &suspend_point);
    // resume immediately if the callback could not be registered
    if (err != hipSuccess) {
      tbb::task::resume(suspend_point);
    }
  });
  HIP_ERROR_CHECK(err);
}
}  // namespace traccc::hip
