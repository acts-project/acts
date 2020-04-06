#pragma once

#include "Acts/Utilities/Platforms/CUDA/CUDAArray.cu"
#include "Acts/Utilities/Platforms/CUDA/CPUMatrix.hpp"

namespace Acts{

template<typename Var_t>
class CUDAMatrix{

public:

  CUDAMatrix()=default;
  CUDAMatrix(size_t nRows, size_t nCols){
    fNRows = nRows;
    fNCols = nCols;
    fSize  = fNRows*fNCols;
    cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t));
  }

  CUDAMatrix(size_t nRows, size_t nCols, CPUMatrix<Var_t>* mat){
    fNRows = nRows;
    fNCols = nCols;
    fSize  = fNRows*fNCols;
    cudaMalloc((Var_t**)&fDevPtr, fNRows*fNCols*sizeof(Var_t));
    CopyH2D(mat->GetEl(0,0),fNRows*fNCols,0);
  }
  
  ~CUDAMatrix(){
    cudaFree(fDevPtr);
  }

  size_t GetNCols(){ return fNCols; }
  size_t GetNRows(){ return fNRows; }
  size_t GetSize() { return fSize; }
  
  Var_t* GetEl(size_t row, size_t col){
    int offset = row+col*fNRows;
    return fDevPtr+offset;
  }

  Var_t* GetHostArray(size_t len, size_t row, size_t col){
    int offset = row+col*fNRows;
    Var_t* hostPtr = new Var_t[len];
    cudaMemcpy(hostPtr, fDevPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);   
    return hostPtr;
  }
  
  CPUArray<Var_t>* GetCPUArray(size_t len, size_t row, size_t col){
    int offset = row+col*fNRows;
    CPUArray<Var_t>* cpuArray = new CPUArray<Var_t>(len);
    cudaMemcpy(cpuArray->Get(), fDevPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);   
    return cpuArray;
  }
  
  void CopyH2D(Var_t* array, size_t len, size_t offset=0){
    cudaMemcpy(fDevPtr+offset, array, len*sizeof(Var_t), cudaMemcpyHostToDevice);
  }

  void CopyH2D(const Var_t* array, size_t len, size_t offset=0){
    cudaMemcpy(fDevPtr+offset, array, len*sizeof(Var_t), cudaMemcpyHostToDevice);
  }

private:
  Var_t* fDevPtr; 
  size_t fNCols;
  size_t fNRows;
  size_t fSize;
};

}

