// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"
#include "Acts/Plugins/Cuda/Utilities/ResultScalar.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace Acts {
namespace Cuda {

template <typename T>
ResultScalar<T>::ResultScalar() : m_array(make_device_array<T>(1)) {}

template <typename T>
typename ResultScalar<T>::pointer ResultScalar<T>::getPtr() const {
  return m_array.get();
}

template <typename T>
typename ResultScalar<T>::Variable_t ResultScalar<T>::get() const {
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

}  // namespace Cuda
}  // namespace Acts

/// Helper macro for instantiating the template code for a given type
#define INST_RSCALAR_FOR_TYPE(TYPE) \
  template class Acts::Cuda::ResultScalar<TYPE>

// Instantiate the templated functions for all primitive types.
INST_RSCALAR_FOR_TYPE(void*);
INST_RSCALAR_FOR_TYPE(char);
INST_RSCALAR_FOR_TYPE(unsigned char);
INST_RSCALAR_FOR_TYPE(short);
INST_RSCALAR_FOR_TYPE(unsigned short);
INST_RSCALAR_FOR_TYPE(int);
INST_RSCALAR_FOR_TYPE(unsigned int);
INST_RSCALAR_FOR_TYPE(long);
INST_RSCALAR_FOR_TYPE(unsigned long);
INST_RSCALAR_FOR_TYPE(long long);
INST_RSCALAR_FOR_TYPE(unsigned long long);
INST_RSCALAR_FOR_TYPE(float);
INST_RSCALAR_FOR_TYPE(double);

// Clean up.
#undef INST_RSCALAR_FOR_TYPE
