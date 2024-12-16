// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

#include <cuda_runtime_api.h>

namespace Acts::detail {

inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    // throw std::runtime_error(ss.str());
    std::cout << ss.str() << std::endl;
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

template <class T>
__global__ void rescaleFeature(std::size_t nbHits, T *data, T scale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  data[i] *= scale;
}

constexpr float g_pi = 3.141592654f;

template <typename T>
__device__ T resetAngle(T angle) {
  if (angle > g_pi) {
    return angle - 2.f * g_pi;
  }
  if (angle < -g_pi) {
    return angle + 2.f * g_pi;
  }
  return angle;
};

template <typename T>
__global__ void makeEdgeFeatures(std::size_t nEdges, const int *srcEdges,
                                 const int *tgtEdges, std::size_t nNodeFeatures,
                                 const T *nodeFeatures, T *edgeFeatures) {
  enum NodeFeatures { r = 0, phi, z, eta };
  constexpr static int nEdgeFeatures = 6;

  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nEdges) {
    return;
  }

  const int src = srcEdges[i];
  const int tgt = tgtEdges[i];

  const T *srcNodeFeatures = nodeFeatures + src * nNodeFeatures;
  const T *tgtNodeFeatures = nodeFeatures + tgt * nNodeFeatures;

  T dr = tgtNodeFeatures[r] - srcNodeFeatures[r];
  T dphi =
      resetAngle(g_pi * (tgtNodeFeatures[phi] - srcNodeFeatures[phi])) / g_pi;
  T dz = tgtNodeFeatures[z] - srcNodeFeatures[z];
  T deta = tgtNodeFeatures[eta] - srcNodeFeatures[eta];
  T phislope = 0.0;
  T rphislope = 0.0;

  if (dr != 0.0) {
    phislope = std::clamp(dphi / dr, -100.f, 100.f);
    T avgR = T{0.5} * (tgtNodeFeatures[r] + srcNodeFeatures[r]);
    rphislope = avgR * phislope;
  }

  T *efPtr = edgeFeatures + i * nEdgeFeatures;
  efPtr[0] = dr;
  efPtr[1] = dphi;
  efPtr[2] = dz;
  efPtr[3] = deta;
  efPtr[4] = phislope;
  efPtr[5] = rphislope;
}

__global__ void setHitId(std::size_t nHits, std::uint64_t *hit_ids) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nHits) {
    return;
  }
  hit_ids[i] = i;
}

__global__ void remapEdges(std::size_t nEdges, int *srcNodes, int *tgtNodes,
                           const std::uint64_t *hit_ids, std::size_t nAllNodes,
                           std::size_t nCompressedNodes) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nEdges) {
    return;
  }

  srcNodes[i] = hit_ids[srcNodes[i]];
  tgtNodes[i] = hit_ids[tgtNodes[i]];
}

template <typename T>
void copyFromDeviceAndPrint(T *data, std::size_t size, std::string_view name) {
  std::vector<T> data_cpu(size);
  CUDA_CHECK(cudaMemcpy(data_cpu.data(), data, size * sizeof(T),
                        cudaMemcpyDeviceToHost));
  std::cout << name << "[" << size << "]: ";
  for (int i = 0; i < size; ++i) {
    std::cout << data_cpu.at(i) << "  ";
  }
  std::cout << std::endl;
}

template <class T>
__global__ void computeXandY(std::size_t nbHits, T *cuda_x, T *cuda_y,
                             const T *cuda_R, const T *cuda_phi) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  double r = cuda_R[i];
  double phi = cuda_phi[i];

  cuda_x[i] = r * std::cos(phi);
  cuda_y[i] = r * std::sin(phi);
}

}  // namespace Acts::detail

#define CUDA_CHECK(ans)                                  \
  do {                                                   \
    Acts::detail::cudaAssert((ans), __FILE__, __LINE__); \
  } while (0)

#define CUDA_PRINTV(ptr, size)             \
  if (logger().level() == Acts::VERBOSE) { \
    std::cout << #ptr << ": ";             \
    cudaDeviceSynchronize();               \
    cudaPrintArray<<<1, 1>>>(ptr, size);   \
    cudaDeviceSynchronize();               \
    std::cout << std::endl;                \
  }
