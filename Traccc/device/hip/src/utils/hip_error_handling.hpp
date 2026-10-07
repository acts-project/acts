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

/// Helper macro used for checking @c hipError_t type return values.
#define TRACCC_HIP_ERROR_CHECK(EXP)                                           \
  do {                                                                        \
    hipError_t errorCode = EXP;                                               \
    if (errorCode != hipSuccess) {                                            \
      traccc::hip::details::throw_error(errorCode, #EXP, __FILE__, __LINE__); \
    }                                                                         \
  } while (false)

namespace traccc::hip::details {

/// Function used to print and throw a user-readable error if something breaks
void throw_error(hipError_t errorCode, const char* expression, const char* file,
                 int line);

}  // namespace traccc::hip::details
