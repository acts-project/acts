// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/device/concepts/barrier.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {

/// A HIP barrier implementation
struct barrier {
  __device__ inline void blockBarrier() const { __syncthreads(); }

  __device__ inline bool blockAnd(bool predicate) const {
    return __syncthreads_and(predicate);
  }

  __device__ inline bool blockOr(bool predicate) const {
    return __syncthreads_or(predicate);
  }

  __device__ inline int blockCount(bool predicate) const {
    return __syncthreads_count(predicate);
  }

};  // struct barrier

// Ensure that the barrier concept is satisfied
static_assert(traccc::device::concepts::barrier<barrier>);

}  // namespace traccc::hip
