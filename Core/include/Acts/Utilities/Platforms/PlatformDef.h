#pragma once
#ifndef PLATFORMDEF
#define PLATFORMDEF

namespace Acts{
  class CPU;
}

#endif

#ifdef ACTS_HAS_CUDA

#include "Acts/Utilities/Platforms/CUDA/CudaScalar.cu"
#include "Acts/Utilities/Platforms/CUDA/CudaVector.cu"
#include "Acts/Utilities/Platforms/CUDA/CudaMatrix.cu"
#include "Acts/Utilities/Platforms/CUDA/CpuMatrix.hpp"

#define WARP_SIZE 32
#define MAX_BLOCK_SIZE 1024

namespace Acts{
  class CUDA;
}

#endif



