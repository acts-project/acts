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
  
template <typename T>
class DeviceMatrix {

public:
  /// The variable type being used
  typedef T Variable_t;

  /// Create a matrix in device memory
  DeviceMatrix(size_t nRows, size_t nCols);

  /// Get the number of rows in the matrix
  size_t nRows() const { return m_nRows; }
  /// Get the number of columns in the matrix
  size_t nCols() const { return m_nCols; }
  /// Get the total size of the matrix
  size_t size() const { return m_nRows * m_nCols; }

  /// Get a specific element of the matrix. (non-const)
  Variable_t& get(size_t row = 0, size_t col = 0);
  /// Get a specific element of the matrix. (const)
  const Variable_t& get(size_t row = 0, size_t col = 0) const;

  /// Get a "pointer into the matrix" (non-const)
  Variable_t* getPtr(size_t row = 0, size_t col = 0);
  /// Get a "pointer into the matrix" (const)
  const Variable_t* getPtr(size_t row = 0, size_t col = 0) const;

  /// Set a specific element of the matrix
  void set(size_t row, size_t col, Variable_t val);

  /// Copy memory from the host.
  void copyFrom(const Variable_t* hostPtr, size_t len, size_t offset);
  /// Copy memory from the host asynchronously.
  void copyFrom(const Variable_t* hostPtr, size_t len, size_t offset,
                cudaStream_t stream);

  /// Reset the matrix to all zeros
  void zeros();
  
private:
  /// Rows in the matrix
  size_t m_nRows;
  /// Culumns in the matrix
  size_t m_nCols;
  /// Smart pointer managing the matrix's memory
  device_array<Variable_t> m_array;

}; // class DeviceMatrix

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "DeviceMatrix.ipp"
