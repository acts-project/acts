#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "Acts/Utilities/Platforms/CUDA/CudaVector.cu"
#include "Acts/Utilities/Platforms/CUDA/CpuMatrix.hpp"
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
  
  Var_t* Get(size_t row, size_t col){
    int offset = row+col*fNRows;
    return fDevPtr+offset;
  }

  void CopyH2D(Var_t* matrix, size_t len, size_t offset=0){
    cudaErrChk( cudaMemcpy(fDevPtr+offset, matrix, len*sizeof(Var_t), cudaMemcpyHostToDevice) );
  }

  void CopyH2D(const Var_t* matrix, size_t len, size_t offset=0){
    cudaErrChk( cudaMemcpy(fDevPtr+offset, matrix, len*sizeof(Var_t), cudaMemcpyHostToDevice) );
  }

private:
  Var_t* fDevPtr; 
  size_t fNCols;
  size_t fNRows;
  size_t fSize;
};

}

