// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "Acts/Plugins/Cuda/Utilities/CpuScalar.hpp"
#include "CudaUtils.cu"

namespace Acts{

template<typename Var_t>
class CudaScalar{

public:
  CudaScalar(){
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t)) );
  }

  CudaScalar(Var_t* scalar){
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t)) );
    cudaErrChk( cudaMemcpy(fDevPtr, scalar, sizeof(Var_t), cudaMemcpyHostToDevice) );
  }

  CudaScalar(const Var_t* scalar){
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t)) );
    cudaErrChk( cudaMemcpy(fDevPtr, scalar, sizeof(Var_t), cudaMemcpyHostToDevice) );
  }
  
  ~CudaScalar(){ 
    cudaErrChk( cudaFree(fDevPtr) );
  }

  Var_t* Get() { return fDevPtr; }

  Var_t GetHost() {
    Var_t* fHostPtr = new Var_t[1];
    cudaErrChk( cudaMemcpy(fHostPtr, fDevPtr, sizeof(Var_t), cudaMemcpyDeviceToHost) );
    return fHostPtr;
  }

  void Zeros() { cudaMemset(fDevPtr,0,sizeof(Var_t)); }
  
  private:
  Var_t* fDevPtr;  
};
}
