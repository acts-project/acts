/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// CUDA include(s).
#include <cuda_runtime_api.h>

/// Helper macro used for checking @c cudaError_t type return values.
#define TRACCC_CUDA_ERROR_CHECK(EXP)                                      \
    do {                                                                  \
        cudaError_t errorCode = EXP;                                      \
        if (errorCode != cudaSuccess) {                                   \
            traccc::cuda_utils::details::throw_error(errorCode, #EXP,     \
                                                     __FILE__, __LINE__); \
        }                                                                 \
    } while (false)

namespace traccc::cuda_utils::details {

/// Function used to print and throw a user-readable error if something breaks
void throw_error(cudaError_t errorCode, const char* expression,
                 const char* file, int line);

}  // namespace traccc::cuda_utils::details
