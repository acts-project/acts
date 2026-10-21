// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// CUDA include(s).
#include <cuda_runtime_api.h>

/// Helper macro used for checking @c cudaError_t type return values.
#define TRACCC_CUDA_ERROR_CHECK(EXP)                                      \
  do {                                                                    \
    cudaError_t errorCode = EXP;                                          \
    if (errorCode != cudaSuccess) {                                       \
      traccc::cuda_utils::details::throw_error(errorCode, #EXP, __FILE__, \
                                               __LINE__);                 \
    }                                                                     \
  } while (false)

namespace traccc::cuda_utils::details {

/// Function used to print and throw a user-readable error if something breaks
void throw_error(cudaError_t errorCode, const char* expression,
                 const char* file, int line);

}  // namespace traccc::cuda_utils::details
