// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(__CUDACC__)

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

/// Number of threads per Warp
#define WARP_SIZE 32

/// Helper macro used for checking @c cudaError_t type return values.
#define DETRAY_CUDA_ERROR_CHECK(ans)       \
  {                                        \
    cudaAssert((ans), __FILE__, __LINE__); \
  }
inline void cudaAssert(cudaError_t code, const char *file, int line,
                       bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "CUDAassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

#endif
