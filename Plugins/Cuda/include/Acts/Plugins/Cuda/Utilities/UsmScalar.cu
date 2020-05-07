// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "CudaUtils.cu"

namespace Acts{

template<typename var_t>
class UsmScalar{

public:
  UsmScalar(){
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged((var_t**)&m_devPtr, sizeof(var_t)) );
    cudaDeviceSynchronize();
  }

  UsmScalar(var_t scalar){
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged((var_t**)&m_devPtr, sizeof(var_t)) );
    cudaDeviceSynchronize();
    m_devPtr[0]=scalar;
  }

  UsmScalar(var_t* scalar){
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged((var_t**)&m_devPtr, sizeof(var_t)) );
    cudaDeviceSynchronize();
    m_devPtr[0]=*scalar;
  }

  UsmScalar(const var_t* scalar){
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged((var_t**)&m_devPtr, sizeof(var_t)) );
    cudaDeviceSynchronize();
    m_devPtr[0]=*scalar;
  }
  
  ~UsmScalar(){
    cudaDeviceSynchronize();
    ACTS_CUDA_ERROR_CHECK( cudaFree(m_devPtr) );
  }

  var_t* get() { return m_devPtr; }
  void set(var_t scalar) { m_devPtr[0]=scalar; }
  
  private:
  var_t* m_devPtr;  
};
}
