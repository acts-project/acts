#pragma once
#ifndef PLATFORMDEF
#define PLATFORMDEF

#include "Acts/Utilities/Platforms/CUDA/CUDAArray.cu"
#include "Acts/Utilities/Platforms/CUDA/CPUArray.hpp"
#include "Acts/Utilities/Platforms/CUDA/CUDAMatrix.cu"
#include "Acts/Utilities/Platforms/CUDA/CPUMatrix.hpp"

// Type definition for each platform

namespace Acts{

class CPU;
class CUDA;

}

#endif
