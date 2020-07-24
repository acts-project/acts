// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"

// CUDA include(s).
#include <cuda.h>

namespace Acts {
namespace Cuda {

template< typename T >
managed_array< T >
make_managed_array( std::size_t size ) {

  // Allocate the memory.
  T* ptr = nullptr;
  if( size != 0 ) {
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged( &ptr, size * sizeof( T ) ) );
  }
  // Create the smart pointer.
  return managed_array< T >( ptr );
}

template< typename T >
device_array< T >
make_device_array( std::size_t size ) {

  // Allocate the memory.
  T* ptr = nullptr;
  if( size != 0 ) {
    ACTS_CUDA_ERROR_CHECK( cudaMalloc( &ptr, size * sizeof( T ) ) );
  }
  // Create the smart pointer.
  return device_array< T >( ptr );
}

template< typename T >
host_array< T >
make_host_array( std::size_t size ) {

  // Allocate the memory.
  T* ptr = nullptr;
  if( size != 0 ) {
    ACTS_CUDA_ERROR_CHECK( cudaMallocHost( &ptr, size * sizeof( T ) ) );
  }
  // Create the smart pointer.
  return host_array< T >( ptr );
}

} // namespace Cuda
} // namespace Acts
