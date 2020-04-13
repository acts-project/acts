#pragma once

namespace Acts{
  class CPU;
}

#ifdef ACTS_HAS_CUDA

#include "Acts/Utilities/Platforms/CUDA/CudaScalar.cu"
#include "Acts/Utilities/Platforms/CUDA/CudaVector.cu"
#include "Acts/Utilities/Platforms/CUDA/CudaMatrix.cu"
#include "Acts/Utilities/Platforms/CUDA/CpuMatrix.hpp"

namespace Acts{
  class CUDA;
}

#endif



