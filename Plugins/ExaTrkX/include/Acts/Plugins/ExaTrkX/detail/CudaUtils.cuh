// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>

#include <cuda_runtime_api.h>

namespace Acts::detail {

inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    throw std::runtime_error(ss.str());
  }
  cudaDeviceSynchronize();
}

inline __global__ void cudaPrintArray(bool *vector, int size) {
  for (int i = 0; i < size; i++)
    printf("%d ", vector[i]);
}
inline __global__ void cudaPrintArray(int *vector, int size) {
  for (int i = 0; i < size; i++)
    printf("%d ", vector[i]);
}
inline __global__ void cudaprintArray(float *vector, int size) {
  for (int i = 0; i < size; ++i)
    printf("%f ", vector[i]);
}
inline __global__ void cudaPrintArray(std::uint64_t *vector, int size) {
  for (int i = 0; i < size; ++i)
    printf("%lu ", vector[i]);
}

}  // namespace Acts::detail

#define CUDA_CHECK(ans)                                  \
  do {                                                   \
    Acts::detail::cudaAssert((ans), __FILE__, __LINE__); \
  } while (0)

#define CUDA_PRINTV(ptr, size)           \
  do {                                   \
    std::cout << #ptr << ": ";           \
    CUDA_CHECK(cudaDeviceSynchronize()); \
    cudaPrintArray<<<1, 1>>>(ptr, size); \
    CUDA_CHECK(cudaGetLastError());      \
    CUDA_CHECK(cudaDeviceSynchronize()); \
    std::cout << std::endl;              \
  } while (0)
