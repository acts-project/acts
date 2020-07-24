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

template <typename T>
HostVector<T>::HostVector(size_t size)
: m_size(size), m_array(make_host_array<T>(size)) {

}

template <typename T>
typename HostVector<T>::Variable_t&
HostVector<T>::get(size_t offset) {

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get()[offset];
}

template <typename T>
const typename HostVector<T>::Variable_t&
HostVector<T>::get(size_t offset) const {

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get()[offset];
}

template <typename T>
typename HostVector<T>::Variable_t*
HostVector<T>::getPtr(size_t offset) {

  // If the vector is empty, return a null pointer.
  if(m_size == 0) {
    return nullptr;
  }

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get() + offset;
}

template <typename T>
const typename HostVector<T>::Variable_t*
HostVector<T>::getPtr(size_t offset) const {

  // If the vector is empty, return a null pointer.
  if(m_size == 0) {
    return nullptr;
  }

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get() + offset;
}

template <typename T>
void HostVector<T>::set(size_t offset, Variable_t val) {

  // Some security check(s).
  assert(offset < m_size);

  // Set the requested element.
  m_array.get()[offset] = val;
  return;
}

template <typename T>
void HostVector<T>::copyFrom(const Variable_t* devPtr, size_t len,
                             size_t offset) {

  // Some security check(s).
  assert(offset + len <= m_size);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_array.get() + offset, devPtr,
                                   len * sizeof(Variable_t),
                                   cudaMemcpyDeviceToHost));
  return;
}

template <typename T>
void HostVector<T>::copyFrom(const Variable_t* devPtr, size_t len,
                             size_t offset, cudaStream_t stream) {

  // Some security check(s).
  assert(offset + len <= m_size);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpyAsync(m_array.get() + offset, devPtr,
                                        len * sizeof(Variable_t),
                                        cudaMemcpyDeviceToHost, stream));
  return;
}

template <typename T>
void HostVector<T>::zeros() {

  memset(m_array.get(), 0, m_size * sizeof(Variable_t));
  return;
}

} // namespace Cuda
} // namespace Acts
