#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

#define WARP_SIZE 32
#define MAX_BLOCK_SIZE 1024

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}


void SetDevice(int devID, std::string& devName){
  cudaDeviceReset();
  
  cudaError_t error;
  cudaDeviceProp deviceProp;
  error = cudaSetDevice(devID);
  error = cudaGetDevice(&devID);
  
  if (error != cudaSuccess){
    printf("cudaGetDevice returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
  }
  error = cudaGetDeviceProperties(&deviceProp, devID);
  
  if (deviceProp.computeMode == cudaComputeModeProhibited){
    fprintf(stderr, "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n");
    exit(EXIT_SUCCESS);
  }
  
  if (error != cudaSuccess)  {
    printf("cudaGetDeviceProperties returned error %s (code %d), line(%d)\n", cudaGetErrorString(error), error, __LINE__);
  }
  else{
    printf("\n GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);
  }
  
  devName = deviceProp.name;
}
