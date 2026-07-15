/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// HIP include(s).
#include <hip/hip_runtime_api.h>

/// Helper macro used for checking @c hipError_t type return values.
#define TRACCC_HIP_ERROR_CHECK(EXP)                                      \
    do {                                                                 \
        hipError_t errorCode = EXP;                                      \
        if (errorCode != hipSuccess) {                                   \
            traccc::hip::details::throw_error(errorCode, #EXP, __FILE__, \
                                              __LINE__);                 \
        }                                                                \
    } while (false)

namespace traccc::hip::details {

/// Function used to print and throw a user-readable error if something breaks
void throw_error(hipError_t errorCode, const char* expression, const char* file,
                 int line);

}  // namespace traccc::hip::details
