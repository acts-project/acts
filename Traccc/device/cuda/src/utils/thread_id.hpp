/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/device/concepts/thread_id.hpp"

namespace traccc::cuda::details {

/// A CUDA thread identifier type
struct thread_id1 {
    __device__ thread_id1() {}

    inline unsigned int __device__ getLocalThreadId() const {
        return threadIdx.x;
    }

    inline unsigned int __device__ getLocalThreadIdX() const {
        return threadIdx.x;
    }

    inline unsigned int __device__ getGlobalThreadId() const {
        return threadIdx.x + blockIdx.x * blockDim.x;
    }

    inline unsigned int __device__ getGlobalThreadIdX() const {
        return threadIdx.x + blockIdx.x * blockDim.x;
    }

    inline unsigned int __device__ getBlockIdX() const { return blockIdx.x; }

    inline unsigned int __device__ getBlockDimX() const { return blockDim.x; }

    inline unsigned int __device__ getGridDimX() const { return gridDim.x; }

};  // struct thread_id1

/// Verify that @c traccc::cuda::details::thread_id1 fulfills the
/// @c traccc::device::concepts::thread_id1 concept.
static_assert(traccc::device::concepts::thread_id1<thread_id1>);

}  // namespace traccc::cuda::details
