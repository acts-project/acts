// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding/Types.hpp"
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"
#include "Acts/Plugins/Cuda/Utilities/HostVector.hpp"
#include "StreamHandlers.cuh"

// CUDA include(s).
#include <cuda_runtime.h>

// System include(s).
#include <cassert>
#include <cstring>

namespace Acts {
namespace Cuda {

template <typename T>
HostVector<T>::HostVector(std::size_t size)
    : m_size(size), m_array(make_host_array<T>(size)) {}

template <typename T>
typename HostVector<T>::element_reference HostVector<T>::get(
    std::size_t offset) {
  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get()[offset];
}

template <typename T>
typename HostVector<T>::element_const_reference HostVector<T>::get(
    std::size_t offset) const {
  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get()[offset];
}

template <typename T>
typename HostVector<T>::pointer HostVector<T>::getPtr(std::size_t offset) {
  // If the vector is empty, return a null pointer.
  if (m_size == 0) {
    return nullptr;
  }

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get() + offset;
}

template <typename T>
typename HostVector<T>::const_pointer HostVector<T>::getPtr(
    std::size_t offset) const {
  // If the vector is empty, return a null pointer.
  if (m_size == 0) {
    return nullptr;
  }

  // Some security check(s).
  assert(offset < m_size);

  // Return the requested element.
  return m_array.get() + offset;
}

template <typename T>
void HostVector<T>::set(std::size_t offset, Variable_t val) {
  // Some security check(s).
  assert(offset < m_size);

  // Set the requested element.
  m_array.get()[offset] = val;
  return;
}

template <typename T>
void HostVector<T>::copyFrom(const_pointer devPtr, std::size_t len,
                             std::size_t offset) {
  // Some security check(s).
  assert(offset + len <= m_size);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_array.get() + offset, devPtr,
                                   len * sizeof(Variable_t),
                                   cudaMemcpyDeviceToHost));
  return;
}

template <typename T>
void HostVector<T>::copyFrom(const_pointer devPtr, std::size_t len,
                             std::size_t offset,
                             const StreamWrapper& streamWrapper) {
  // Some security check(s).
  assert(offset + len <= m_size);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(
      cudaMemcpyAsync(m_array.get() + offset, devPtr, len * sizeof(Variable_t),
                      cudaMemcpyDeviceToHost, getStreamFrom(streamWrapper)));
  return;
}

template <typename T>
void HostVector<T>::zeros() {
  memset(m_array.get(), 0, m_size * sizeof(Variable_t));
  return;
}

}  // namespace Cuda
}  // namespace Acts

/// Helper macro for instantiating the template code for a given type
#define INST_HVECTOR_FOR_TYPE(TYPE) template class Acts::Cuda::HostVector<TYPE>

// Instantiate the templated functions for all primitive types.
INST_HVECTOR_FOR_TYPE(void*);
INST_HVECTOR_FOR_TYPE(char);
INST_HVECTOR_FOR_TYPE(unsigned char);
INST_HVECTOR_FOR_TYPE(short);
INST_HVECTOR_FOR_TYPE(unsigned short);
INST_HVECTOR_FOR_TYPE(int);
INST_HVECTOR_FOR_TYPE(unsigned int);
INST_HVECTOR_FOR_TYPE(long);
INST_HVECTOR_FOR_TYPE(unsigned long);
INST_HVECTOR_FOR_TYPE(long long);
INST_HVECTOR_FOR_TYPE(unsigned long long);
INST_HVECTOR_FOR_TYPE(float);
INST_HVECTOR_FOR_TYPE(double);

// Instantiate them for any necessary custom type(s) as well.
INST_HVECTOR_FOR_TYPE(Acts::Cuda::details::Triplet);

// Clean up.
#undef INST_HVECTOR_FOR_TYPE
