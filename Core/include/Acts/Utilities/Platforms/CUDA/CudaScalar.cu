#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
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
  
  private:
  Var_t* fDevPtr;  
};
}
