// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/Arrays.cuh"

// CUDA include(s).
#include <cuda.h>

namespace Acts {
namespace Cuda {

/// Vector holding data in device memory
template <typename T>
class DeviceVector {

public:
  /// The variable type being used
  typedef T Variable_t;

  /// Create a vector in the/a device's memory
  DeviceVector(size_t size);

  /// Get the size of the vector
  size_t size() const { return m_size; }

  /// Get a specific element of the vector (non-const)
  Variable_t& get(size_t offset = 0);
  /// Get a specific element of the vector (const)
  const Variable_t& get(size_t offset = 0) const;

  /// Get a "pointer into the vector" (non-const)
  Variable_t* getPtr(size_t offset = 0);
  /// Get a "pointer into the vector" (const)
  const Variable_t* getPtr(size_t offset = 0) const;

  /// Set a specific element of the vector
  void set(size_t offset, Variable_t val);

  /// Copy memory from the host.
  void copyFrom(const Variable_t* hostPtr, size_t len, size_t offset);
  /// Copy memory from the host asynchronously.
  void copyFrom(const Variable_t* hostPtr, size_t len, size_t offset,
                cudaStream_t stream);

  /// Reset the vector to all zeros
  void zeros();

private:
  /// The size of the vector
  size_t m_size;
  /// Smart pointer managing the vector's memory
  device_array< Variable_t > m_array;

}; // class DeviceVector

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "DeviceVector.ipp"
