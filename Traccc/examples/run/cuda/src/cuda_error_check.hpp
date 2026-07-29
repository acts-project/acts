/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <exception>

// Helper macro for checking the return value of CUDA function calls
#define CUDA_ERROR_CHECK(EXP)                                            \
  do {                                                                   \
    const cudaError_t errorCode = EXP;                                   \
    if (errorCode != cudaSuccess) {                                      \
      throw std::runtime_error(std::string("Failed to run " #EXP " (") + \
                               cudaGetErrorString(errorCode) + ")");     \
    }                                                                    \
  } while (false)
