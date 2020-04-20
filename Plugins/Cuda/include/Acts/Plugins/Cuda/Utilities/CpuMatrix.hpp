// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaMatrix.cu"

// column-major style Matrix Definition

namespace Acts{

template<typename Var_t>
class CudaMatrix;
  
template<typename Var_t>
class CpuMatrix{
  
public:

  CpuMatrix() = default;
  CpuMatrix(size_t nRows, size_t nCols, bool pinned=0){
    SetSize(nRows,nCols);
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[fSize];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, fSize*sizeof(Var_t));
    }
  }

  CpuMatrix(size_t nRows, size_t nCols, CudaMatrix<Var_t>* cuMat, bool pinned=0){
    SetSize(nRows,nCols);
    fPinned = pinned;
    if (pinned == 0){
      fHostPtr = new Var_t[fSize];
    }
    else if (pinned == 1){
      cudaMallocHost(&fHostPtr, fNRows*fNCols*sizeof(Var_t));
    }
    cudaMemcpy(fHostPtr, cuMat->Get(0,0), fSize*sizeof(Var_t), cudaMemcpyDeviceToHost);   
  }
  
  ~CpuMatrix(){
    if (!fPinned){
      delete fHostPtr;
    }
    else if (fPinned){
      cudaFreeHost(fHostPtr);
    }
  }

  void SetSize(size_t row, size_t col){
    fNRows = row;
    fNCols = col;
    fSize  = fNRows*fNCols; 
  }
  
  size_t GetNCols(){ return fNCols; }
  size_t GetNRows(){ return fNRows; }
  size_t GetSize() { return fSize; }
  
  Var_t* Get(size_t row=0, size_t col=0){
    size_t offset=row+col*fNRows;
    return fHostPtr+offset;
  }
  
  void Set(size_t row, size_t col, Var_t val){
    size_t offset=row+col*fNRows;
    fHostPtr[offset] = val;
  }
  
  Var_t* GetColumn(size_t col){
    return fHostPtr+col*fNRows;
  }
  Var_t* GetRow(size_t row){    
    Var_t* ret = new Var_t[fNCols];
    for(size_t i_c=0; i_c<fNCols; i_c++) ret[i_c] = fHostPtr[row+fNRows*i_c];
    return ret;    
  }

  void SetRow(size_t row, Var_t* input){
    for(size_t i_c=0; i_c<fNCols; i_c++){
      fHostPtr[row+fNRows*i_c]=input[i_c];
    }
  }

  void SetColumn(size_t col, Var_t* input){
    fHostPtr[col*fNRows] = input[0];
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
  size_t fNCols;
  size_t fNRows;
  size_t fSize;
  bool   fPinned;
};
  
}
