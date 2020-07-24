// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA include(s).
#include <cuda_runtime.h>

/// Helper macro used in the CUDA plugin for checking @c cudaError_t type return
/// values.
#define ACTS_CUDA_ERROR_CHECK( EXP )                                           \
  do {                                                                         \
    cudaError_t errorCode = EXP;                                               \
    if( errorCode != cudaSuccess ) {                                           \
      Acts::Cuda::details::throwError( errorCode, #EXP, __FILE__, __LINE__ );  \
    }                                                                          \
  } while( false )

namespace Acts {
namespace Cuda {
namespace details {

  /// Function used to print and throw a user-readable error if something breaks
  void throwError( cudaError_t errorCode, const char* expression,
                   const char* file, int line );

} // namespace details
} // namespace Cuda
} // namespace Acts
