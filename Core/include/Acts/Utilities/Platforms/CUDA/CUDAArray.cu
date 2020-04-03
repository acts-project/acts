#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "Acts/Utilities/Platforms/CUDA/CPUArray.hpp"

namespace Acts{

template<typename Var_t>
class CUDAArray{

public:
  
  CUDAArray(size_t size){ 
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
  }

  CUDAArray(size_t size, Var_t* buffer, size_t len, size_t offset=0){ 
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
    CopyH2D(buffer, len, offset);
  }

  CUDAArray(size_t size, const Var_t* buffer, size_t len, size_t offset=0){ 
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
    CopyH2D(buffer, len, offset);
  }
  
  ~CUDAArray(){ 
    cudaFree(fDevPtr); 
  }

  size_t GetSize(){return fSize;}
  
  Var_t* Get(size_t offset=0) { return fDevPtr+offset; }
  
  Var_t* GetHostArray(size_t len, size_t offset=0) const {
    Var_t* hostArray = new Var_t[len];
    cudaMemcpy(hostArray, fDevPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);   
    return hostArray;
  }

  CPUArray<Var_t>* GetCPUArray(size_t len, size_t offset=0) const {
    CPUArray<Var_t>* cpuArray = new CPUArray<Var_t>(len);
    cudaMemcpy(cpuArray->Get(), fDevPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);   
    return cpuArray;
  }
  
  //Var_t& operator[](std::size_t idx)       { return fDevPtr[idx]; }  // Need to test
  //const Var_t& operator[](std::size_t idx) const { return fDevPtr[idx]; }  // Need to test
  
  void CopyH2D(Var_t* array, size_t len, size_t offset=0){
    cudaMemcpy(fDevPtr+offset, array, len*sizeof(Var_t), cudaMemcpyHostToDevice);
  }

  void CopyH2D(const Var_t* array, size_t len, size_t offset=0){
    cudaMemcpy(fDevPtr+offset, array, len*sizeof(Var_t), cudaMemcpyHostToDevice);
  }
    
private:
  Var_t* fDevPtr; 
  size_t fSize;
};
}
