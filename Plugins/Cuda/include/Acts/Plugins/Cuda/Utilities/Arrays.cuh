// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <memory>

namespace Acts {
namespace Cuda {

/// Namespace holding some implementation detail types that should not be used
/// directly.
namespace details {

/// Class performing the deletion of a CUDA managed memory array
class ManagedArrayDeleter {

public:
  /// Operator performing the deletion of the memory
  void operator()( void* ptr );

}; // class ManagedArrayDeleter

/// Class performing the deletion of a CUDA device memory array
class DeviceArrayDeleter {

public:
  /// Operator performing the deletion of the memory
  void operator()( void* ptr );

}; // class DeviceArrayDeleter

/// Class performing the deletion of pinned host memory
class HostArrayDeleter {

public:
  /// Operator performing the deletion of the memory
  void operator()( void* ptr );

}; // class HostArrayDeleter

} // namespace details

/// Convenience type for handling primitive variable arrays in CUDA managed
/// memory
template< typename T >
using managed_array = std::unique_ptr< T, details::ManagedArrayDeleter >;

/// Function creating a primitive array in CUDA managed memory
template< typename T >
managed_array< T > make_managed_array( std::size_t size );

/// Convenience type for using primitive variable arrays on a CUDA device
template< typename T >
using device_array = std::unique_ptr< T, details::DeviceArrayDeleter >;

/// Function creating a primitive array in CUDA device memory
template< typename T >
device_array< T > make_device_array( std::size_t size );

/// Convenience type for using primitive variable arrays on the host
template< typename T >
using host_array = std::unique_ptr< T, details::HostArrayDeleter >;

/// Function creating a primitive array in the host's memory
template< typename T >
host_array< T > make_host_array( std::size_t size );

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "Arrays.ipp"
