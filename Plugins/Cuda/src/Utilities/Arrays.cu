// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA include(s).
#include <cuda.h>

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/Arrays.cuh"
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"

namespace Acts {
namespace Cuda {
namespace details {

  void ManagedArrayDeleter::operator()( void* ptr ) {

    // Ignore null-pointers.
    if( ptr == nullptr ) {
      return;
    }

    // Free the managed memory.
    ACTS_CUDA_ERROR_CHECK( cudaFree( ptr ) );
    return;
  }

  void DeviceArrayDeleter::operator()( void* ptr ) {

    // Ignore null-pointers.
    if( ptr == nullptr ) {
      return;
    }

    // Free the device memory.
    ACTS_CUDA_ERROR_CHECK( cudaFree( ptr ) );
    return;
  }

  void HostArrayDeleter::operator()( void* ptr ) {

    // Ignore null-pointers.
    if( ptr == nullptr ) {
      return;
    }

    // Free the pinned host memory.
    ACTS_CUDA_ERROR_CHECK( cudaFreeHost( ptr ) );
    return;
  }

} // namespace details
} // namespace Cuda
} // namespace Acts
