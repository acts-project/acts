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

namespace Acts {
namespace Cuda {

/// Helper class for handling return scalar values from CUDA kernels
///
/// Note that this is only needed for return values. Input scalar values can be
/// passed simply as function arguments to CUDA kernels.
///
template <typename T>
class ResultScalar {

public:
  /// The variable type being used
  typedef T Variable_t;

  /// Default constructor
  ResultScalar();

  /// Get the "device pointer" that the kernel can write to
  Variable_t* getPtr() const;

  /// Get the value set by the device
  Variable_t get() const;
  /// Convenience operator for using the result value
  operator Variable_t() const;

private:
  /// Variable managing the memory on the device
  device_array<Variable_t> m_array;

}; // class ResultScalar

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "ResultScalar.ipp"
