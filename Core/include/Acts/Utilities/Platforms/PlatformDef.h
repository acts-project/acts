#pragma once
#ifndef PLATFORMDEF
#define PLATFORMDEF

namespace Acts{
  class CPU;
}

#endif

#ifdef ACTS_HAS_CUDA

#include "Acts/Utilities/Platforms/CUDA/CUDAArray.cu"
#include "Acts/Utilities/Platforms/CUDA/CPUArray.hpp"
#include "Acts/Utilities/Platforms/CUDA/CUDAMatrix.cu"
#include "Acts/Utilities/Platforms/CUDA/CPUMatrix.hpp"

#define WARP_SIZE 32
#define MAX_BLOCK_SIZE 1024

namespace Acts{
  class CUDA;
}

#endif
