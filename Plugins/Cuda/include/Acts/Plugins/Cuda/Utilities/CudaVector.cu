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
#include "Acts/Plugins/Cuda/Utilities/CpuVector.hpp"
#include "CudaUtils.cu"

namespace Acts{

template<typename Var_t>
class CudaVector{

public:
  
  CudaVector(size_t size){ 
    fSize = size;
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t)) );
  }

  CudaVector(size_t size, Var_t* vector){
    fSize = size;
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t)) );
    CopyH2D(vector, fSize, 0);
  }
     
  CudaVector(size_t size, Var_t* vector, size_t len, size_t offset){ 
    fSize = size;
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t)) );
    CopyH2D(vector, len, offset);
  }
  
  ~CudaVector(){ 
    cudaFree(fDevPtr); 
  }

  size_t GetSize(){return fSize;}
  
  Var_t* Get(size_t offset=0) { return fDevPtr+offset; }

  Var_t* GetHost() {
    Var_t* fHostPtr = new Var_t[fSize];
    cudaErrChk( cudaMemcpy(fHostPtr, fDevPtr, fSize*sizeof(Var_t), cudaMemcpyDeviceToHost) );
    return fHostPtr;
  }

  void CopyH2D(Var_t* vector, size_t len, size_t offset){
    cudaErrChk( cudaMemcpy(fDevPtr+offset, vector, len*sizeof(Var_t), cudaMemcpyHostToDevice) );
  }
  void CopyH2D(Var_t* vector, size_t len, size_t offset, cudaStream_t* stream){
    cudaErrChk( cudaMemcpyAsync(fDevPtr+offset, vector, len*sizeof(Var_t), cudaMemcpyHostToDevice, *stream) );
  }

  void Zeros(){
    cudaErrChk( cudaMemset(fDevPtr, 0, fSize*sizeof(Var_t)) );
  }
  
private:
  Var_t* fDevPtr; 
  size_t fSize;
};
}
