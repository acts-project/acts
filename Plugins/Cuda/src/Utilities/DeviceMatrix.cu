// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding/Types.hpp"
#include "Acts/Plugins/Cuda/Utilities/DeviceMatrix.hpp"
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"
#include "StreamHandlers.cuh"

// CUDA include(s).
#include <cuda_runtime.h>

// System include(s).
#include <cassert>
#include <cstring>

namespace Acts {
namespace Cuda {

template <typename T>
DeviceMatrix<T>::DeviceMatrix(std::size_t nRows, std::size_t nCols)
    : m_nRows(nRows),
      m_nCols(nCols),
      m_array(make_device_array<T>(nRows * nCols)) {}

template <typename T>
typename DeviceMatrix<T>::element_reference DeviceMatrix<T>::get(
    std::size_t row, std::size_t col) {
  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template <typename T>
typename DeviceMatrix<T>::element_const_reference DeviceMatrix<T>::get(
    std::size_t row, std::size_t col) const {
  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template <typename T>
typename DeviceMatrix<T>::pointer DeviceMatrix<T>::getPtr(std::size_t row,
                                                          std::size_t col) {
  // If the matrix is empty, just return a null pointer.
  if ((m_nRows == 0) || (m_nCols == 0)) {
    return nullptr;
  }

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template <typename T>
typename DeviceMatrix<T>::const_pointer DeviceMatrix<T>::getPtr(
    std::size_t row, std::size_t col) const {
  // If the matrix is empty, just return a null pointer.
  if ((m_nRows == 0) || (m_nCols == 0)) {
    return nullptr;
  }

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template <typename T>
void DeviceMatrix<T>::set(std::size_t row, std::size_t col, Variable_t val) {
  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Set the requested element.
  m_array.get()[row + col * m_nRows] = val;
  return;
}

template <typename T>
void DeviceMatrix<T>::copyFrom(const_pointer hostPtr, std::size_t len,
                               std::size_t offset) {
  // Some security check(s).
  assert(offset + len <= m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_array.get() + offset, hostPtr,
                                   len * sizeof(Variable_t),
                                   cudaMemcpyHostToDevice));
  return;
}

template <typename T>
void DeviceMatrix<T>::copyFrom(const_pointer hostPtr, std::size_t len,
                               std::size_t offset,
                               const StreamWrapper& streamWrapper) {
  // Some security check(s).
  assert(offset + len <= m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(
      cudaMemcpyAsync(m_array.get() + offset, hostPtr, len * sizeof(Variable_t),
                      cudaMemcpyHostToDevice, getStreamFrom(streamWrapper)));
  return;
}

template <typename T>
void DeviceMatrix<T>::zeros() {
  memset(m_array.get(), 0, m_nRows * m_nCols * sizeof(Variable_t));
  return;
}

}  // namespace Cuda
}  // namespace Acts

/// Helper macro for instantiating the template code for a given type
#define INST_DMATRIX_FOR_TYPE(TYPE) \
  template class Acts::Cuda::DeviceMatrix<TYPE>

// Instantiate the templated functions for all primitive types.
INST_DMATRIX_FOR_TYPE(void*);
INST_DMATRIX_FOR_TYPE(char);
INST_DMATRIX_FOR_TYPE(unsigned char);
INST_DMATRIX_FOR_TYPE(short);
INST_DMATRIX_FOR_TYPE(unsigned short);
INST_DMATRIX_FOR_TYPE(int);
INST_DMATRIX_FOR_TYPE(unsigned int);
INST_DMATRIX_FOR_TYPE(long);
INST_DMATRIX_FOR_TYPE(unsigned long);
INST_DMATRIX_FOR_TYPE(long long);
INST_DMATRIX_FOR_TYPE(unsigned long long);
INST_DMATRIX_FOR_TYPE(float);
INST_DMATRIX_FOR_TYPE(double);

// Instantiate them for any necessary custom type(s) as well.
INST_DMATRIX_FOR_TYPE(Acts::Cuda::details::Triplet);

// Clean up.
#undef INST_DMATRIX_FOR_TYPE
