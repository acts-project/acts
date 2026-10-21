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
