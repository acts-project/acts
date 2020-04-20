// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaVector.cu"

namespace Acts{

template<typename Var_t>
class CudaVector;
  
template<typename Var_t>
class CpuVector{
  
public:

  CpuVector() = default;
  CpuVector(size_t size, bool pinned=0){
    fSize = size;
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[fSize];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, fSize*sizeof(Var_t));
    }
  }

  CpuVector(size_t size, CudaVector<Var_t>* cuVec, bool pinned=0){
    fSize = size;
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[fSize];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, fSize*sizeof(Var_t));
    }
    cudaMemcpy(fHostPtr, cuVec->Get(), fSize*sizeof(Var_t), cudaMemcpyDeviceToHost);   
  }
  
  ~CpuVector(){
    if (!fPinned){
      delete fHostPtr;
    }
    else if (fPinned){
      cudaFreeHost(fHostPtr);
    }
  }

  size_t GetSize() { return fSize; }
  
  Var_t* Get(size_t offset=0){
    return fHostPtr+offset;
  }
  
  void Set(size_t offset, Var_t val){
    fHostPtr[offset] = val;
  }
  
  void CopyD2H(Var_t* devPtr, size_t len, size_t offset){
    cudaMemcpy(fHostPtr+offset, devPtr, len*sizeof(Var_t), cudaMemcpyDeviceToHost);
  }

  void CopyD2H(Var_t* devPtr, size_t len, size_t offset, cudaStream_t* stream){
    cudaMemcpyAsync(fHostPtr+offset, devPtr, len*sizeof(Var_t), cudaMemcpyDeviceToHost, *stream);
  }

  void Zeros(){
    memset(fHostPtr, 0, fSize*sizeof(Var_t));
  }
  
private:
  Var_t* fHostPtr;
  size_t fSize;
  bool   fPinned;
};
  
}

