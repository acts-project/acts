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
#include "Acts/Plugins/Cuda/Utilities/CpuMatrix.hpp"
#include "CudaUtils.cu"

namespace Acts{

template<typename Var_t>
class CudaMatrix{

public:

  CudaMatrix()=default;
  CudaMatrix(size_t nRows, size_t nCols){
    SetSize(nRows,nCols);
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t)) );
  }

  CudaMatrix(size_t nRows, size_t nCols, Var_t* mat){
    SetSize(nRows,nCols);
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t)) );
    CopyH2D(mat, fSize, 0);
  }
  
  CudaMatrix(size_t nRows, size_t nCols, CpuMatrix<Var_t>* mat){
    SetSize(nRows,nCols);
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t)) );
    CopyH2D(mat->Get(0,0), fSize, 0);
  }

  CudaMatrix(size_t nRows, size_t nCols, Var_t* mat, size_t len, size_t offset){
    SetSize(nRows,nCols);
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t)) );
    CopyH2D(mat, len, offset);
  }
  
  CudaMatrix(size_t nRows, size_t nCols, CpuMatrix<Var_t>* mat, size_t len, size_t offset){
    SetSize(nRows,nCols);
    cudaErrChk( cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t)) );
    CopyH2D(mat->Get(0,0),len,offset);
  }
  
  ~CudaMatrix(){
    cudaFree(fDevPtr);
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
    int offset = row+col*fNRows;
    return fDevPtr+offset;
  }

  void CopyH2D(Var_t* matrix, size_t len, size_t offset=0){
    cudaErrChk( cudaMemcpy(fDevPtr+offset, matrix, len*sizeof(Var_t), cudaMemcpyHostToDevice) );
  }

  void CopyH2D(const Var_t* matrix, size_t len, size_t offset=0){
    cudaErrChk( cudaMemcpy(fDevPtr+offset, matrix, len*sizeof(Var_t), cudaMemcpyHostToDevice) );
  }
  
  void Zeros(){
    cudaErrChk( cudaMemset(fDevPtr, 0, fSize*sizeof(Var_t)) );
  }
  
private:
  Var_t* fDevPtr; 
  size_t fNCols;
  size_t fNRows;
  size_t fSize;
};

}

