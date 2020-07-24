// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/Arrays.hpp"
#include "Acts/Plugins/Cuda/Utilities/StreamWrapper.hpp"

namespace Acts {
namespace Cuda {

/// Vector holding data in host-pinned memory
template <typename T>
class HostVector {

public:
  /// The variable type being used
  typedef T Variable_t;

  /// Create a vector in host memory
  HostVector(std::size_t size);

  /// Get the size of the vector
  std::size_t size() const { return m_size; }

  /// Get a specific element of the vector (non-const)
  Variable_t& get(std::size_t offset = 0);
  /// Get a specific element of the vector (const)
  const Variable_t& get(std::size_t offset = 0) const;

  /// Get a "pointer into the vector" (non-const)
  Variable_t* getPtr(std::size_t offset = 0);
  /// Get a "pointer into the vector" (const)
  const Variable_t* getPtr(std::size_t offset = 0) const;

  /// Set a specific element of the vector
  void set(std::size_t offset, Variable_t val);

  /// Copy memory from a/the device.
  void copyFrom(const Variable_t* devPtr, std::size_t len, std::size_t offset);
  /// Copy memory from a/the device asynchronously.
  void copyFrom(const Variable_t* devPtr, std::size_t len, std::size_t offset,
                const StreamWrapper& streamWrapper);

  /// Reset the vector to all zeros
  void zeros();

private:
  /// The size of the vector
  std::size_t m_size;
  /// Smart pointer managing the vector's memory
  host_array< Variable_t > m_array;

}; // class HostVector

} // namespace Cuda
} // namespace Acts
