// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaScalar.cu"

namespace Acts{

template<typename Var_t>
class CudaScalar;
  
template<typename Var_t>
class CpuScalar{
  
public:

  CpuScalar() = default;
  CpuScalar(bool pinned=0){
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[1];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, sizeof(Var_t));
    }
  }

  CpuScalar(CudaScalar<Var_t>* cuScalar, bool pinned=0){
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[1];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, sizeof(Var_t));
    }
    cudaMemcpy(fHostPtr, cuScalar->Get(), sizeof(Var_t), cudaMemcpyDeviceToHost);   
  }
  
  ~CpuScalar(){
    if (!fPinned){
      delete fHostPtr;
    }
    else if (fPinned){
      cudaFreeHost(fHostPtr);
    }
  }
  
  Var_t* Get(){
    return fHostPtr;
  }
  
  void Set(Var_t val){
    fHostPtr[0] = val;
  }
    
private:
  Var_t* fHostPtr;
  size_t fSize;
  bool   fPinned;
};
  
}

