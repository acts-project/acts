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

// System include(s).
#include <cassert>
#include <cstring>

namespace Acts {
namespace Cuda {

template< typename T >
HostMatrix< T >::HostMatrix( size_t nRows, size_t nCols )
: m_nRows( nRows ), m_nCols( nCols ),
  m_array( make_host_array< T >( nRows * nCols ) ) {

}

template< typename T >
typename HostMatrix< T >::Variable_t&
HostMatrix< T >::get(size_t row, size_t col) {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template< typename T >
const typename HostMatrix< T >::Variable_t&
HostMatrix< T >::get(size_t row, size_t col) const {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template< typename T >
typename HostMatrix< T >::Variable_t*
HostMatrix< T >::getPtr(size_t row, size_t col) {

  // If the matrix is empty, just return a null pointer.
  if((m_nRows == 0) || (m_nCols == 0)) {
    return nullptr;
  }

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template< typename T >
const typename HostMatrix< T >::Variable_t*
HostMatrix< T >::getPtr(size_t row, size_t col) const {

  // If the matrix is empty, just return a null pointer.
  if((m_nRows == 0) || (m_nCols == 0)) {
    return nullptr;
  }

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template< typename T >
void HostMatrix< T >::set(size_t row, size_t col, Variable_t val) {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Set the requested element.
  m_array.get()[row + col * m_nRows] = val;
  return;
}

template< typename T >
void HostMatrix< T >::copyFrom(const Variable_t* devPtr, size_t len,
                               size_t offset) {

  // Some security check(s).
  assert(offset + len <= m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_array.get() + offset, devPtr,
                                   len * sizeof(Variable_t),
                                   cudaMemcpyDeviceToHost));
  return;
}

template< typename T >
void HostMatrix< T >::copyFrom(const Variable_t* devPtr, size_t len,
                               size_t offset, cudaStream_t stream) {

  // Some security check(s).
  assert(offset + len <= m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpyAsync(m_array.get() + offset, devPtr,
                                        len * sizeof(Variable_t),
                                        cudaMemcpyDeviceToHost, stream));
  return;
}

template< typename T >
void HostMatrix< T >::zeros() {

  memset(m_array.get(), 0, m_nRows * m_nCols * sizeof(Variable_t));
  return;
}

} // namespace Cuda
} // namespace Acts
