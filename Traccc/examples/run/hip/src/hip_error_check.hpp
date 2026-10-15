// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// HIP include(s).
#include <hip/hip_runtime_api.h>

// System include(s).
#include <stdexcept>
#include <string>

// Helper macro for checking the return value of HIP function calls
#define HIP_ERROR_CHECK(EXP)                                             \
  do {                                                                   \
    const hipError_t errorCode = EXP;                                    \
    if (errorCode != hipSuccess) {                                       \
      throw std::runtime_error(std::string("Failed to run " #EXP " (") + \
                               hipGetErrorString(errorCode) + ")");      \
    }                                                                    \
  } while (false)
