#ifndef CPUARRAY
#define CPUARRAY

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "Acts/Utilities/Platforms/CUDA/CUDAArray.cu"

namespace Acts{

template<typename Var_t>
class CUDAArray;
  
template<typename Var_t>
class CPUArray{

public:

  CPUArray()=default;
  
  CPUArray(size_t size){ 
    fSize = size;
    cudaMallocHost(&fHostPtr, fSize*sizeof(Var_t));
  }

  CPUArray(size_t size, CUDAArray<Var_t>* cuBuf){ 
    fSize = size;
    fHostPtr = (cuBuf->GetCPUArray(fSize,0))->Get();
  }
  ~CPUArray(){ cudaFreeHost(fHostPtr); }

  Var_t* Get(size_t offset=0){ return fHostPtr+offset; }

  void CopyD2H(Var_t* devPtr, size_t len, size_t offset=0){
    cudaMemcpy(fHostPtr, devPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);
  }

   	Var_t& operator[](std::size_t idx)       { return fHostPtr[idx]; }
  const Var_t& operator[](std::size_t idx) const { return fHostPtr[idx]; }

  
private:
  Var_t* fHostPtr; 
  size_t fSize;   
};

}

#endif
