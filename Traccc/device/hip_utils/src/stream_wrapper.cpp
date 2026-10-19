// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/hip/utils/stream_wrapper.hpp"

#include "hip_error_handling.hpp"

// HIP include(s).
#include <hip/hip_runtime_api.h>

namespace traccc::hip {

stream_wrapper::stream_wrapper(void* stream) : m_stream{stream} {}

int stream_wrapper::device() const {
  // The device ID.
  int device = -1;

#ifdef TRACCC_HAVE_HIP_STREAM_GET_DEVICE
  // Somewhere around ROCm 6.0 this became available.
  TRACCC_HIP_ERROR_CHECK(
      hipStreamGetDevice(static_cast<hipStream_t>(m_stream), &device));

#else
  // If the HIP version is too old to support hipStreamGetDevice, we return
  // the current device instead. This is not ideal, but should be good enough.
  TRACCC_HIP_ERROR_CHECK(hipGetDevice(&device));
#endif

  // Return the device ID.
  return device;
}

void* stream_wrapper::hipStream() const {
  return m_stream;
}

void stream_wrapper::synchronize() const {
  TRACCC_HIP_ERROR_CHECK(
      hipStreamSynchronize(static_cast<hipStream_t>(m_stream)));
}

}  // namespace traccc::hip
