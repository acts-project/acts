#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"

namespace Acts{

template<typename Var_t>
class CudaScalar{

public:
  CudaScalar(){
    cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t));  
  }

  CudaScalar(Var_t* scalar){
    cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t));
    cudaMemcpy(fDevPtr, scalar, sizeof(Var_t), cudaMemcpyHostToDevice);
  }

  CudaScalar(const Var_t* scalar){
    cudaMalloc((Var_t**)&fDevPtr, sizeof(Var_t));
    cudaMemcpy(fDevPtr, scalar, sizeof(Var_t), cudaMemcpyHostToDevice);
  }
  
  ~CudaScalar(){ 
    cudaFree(fDevPtr); 
  }

  Var_t* Get() { return fDevPtr; }
  
  private:
  Var_t* fDevPtr;  
};
}
