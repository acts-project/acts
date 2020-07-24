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

template <typename T>
ResultScalar<T>::ResultScalar()
: m_array(make_device_array<T>(1)) {

}

template <typename T>
typename ResultScalar<T>::Variable_t*
ResultScalar<T>::getPtr() const {

  return m_array.get();
}

template <typename T>
typename ResultScalar<T>::Variable_t
ResultScalar<T>::get() const {

  // Retrieve the value of the variable.
  Variable_t result;
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(&result, m_array.get(), sizeof(Variable_t),
                        cudaMemcpyDeviceToHost));
  // Return it.
  return result;
}

template <typename T>
ResultScalar<T>::operator typename ResultScalar<T>::Variable_t() const {

  // Rely on the get() function...
  return get();
}

} // namespace Cuda
} // namespace Acts
