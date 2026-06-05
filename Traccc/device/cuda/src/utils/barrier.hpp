/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::cuda {

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
};

}  // namespace traccc::cuda
