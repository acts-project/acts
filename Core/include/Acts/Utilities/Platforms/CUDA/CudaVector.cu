#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"

namespace Acts{

template<typename Var_t>
class CudaVector{

public:
  
  CudaVector(size_t size){ 
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
  }

  CudaVector(size_t size, Var_t* vector){
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
    CopyH2D(vector, fSize, 0);
  }
     
  CudaVector(size_t size, Var_t* vector, size_t len, size_t offset){ 
    fSize = size;
    cudaMalloc((Var_t**)&fDevPtr, fSize*sizeof(Var_t));
    CopyH2D(vector, len, offset);
  }
  
  ~CudaVector(){ 
    cudaFree(fDevPtr); 
  }

  size_t GetSize(){return fSize;}
  
  Var_t* Get(size_t offset=0) { return fDevPtr+offset; }

  Var_t* GetHost() {
    Var_t* fHostPtr = new Var_t[fSize];
    cudaMemcpy(fHostPtr, fDevPtr, fSize*sizeof(Var_t), cudaMemcpyDeviceToHost);
    return fHostPtr;
  }

  void CopyH2D(Var_t* vector, size_t len, size_t offset){
    cudaMemcpy(fDevPtr+offset, vector, len*sizeof(Var_t), cudaMemcpyHostToDevice);
  }
  void CopyH2D(Var_t* vector, size_t len, size_t offset, cudaStream_t* stream){
    cudaMemcpyAsync(fDevPtr+offset, vector, len*sizeof(Var_t), cudaMemcpyHostToDevice, *stream);
  }
      
private:
  Var_t* fDevPtr; 
  size_t fSize;
};
}
