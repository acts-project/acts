// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <cuda_runtime_api.h>

namespace ActsPlugins::detail {

constexpr float g_pi = 3.14159265358979323846f;

template <typename T>
__device__ T resetAngle(T angle) {
  if (angle > g_pi) {
    return angle - 2.f * g_pi;
  }
  if (angle < -g_pi) {
    return angle + 2.f * g_pi;
  }
  assert(angle >= -g_pi && angle < g_pi);
  return angle;
};

constexpr int g_nEdgeFeatures = 6;

template <typename T>
__device__ void computeEdgeFeatures(int src, int tgt, std::size_t nNodeFeatures,
                                    const T *nodeFeatures, T *efPtr) {
  enum NodeFeatures { r = 0, phi, z, eta };

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

  efPtr[0] = dr;
  efPtr[1] = dphi;
  efPtr[2] = dz;
  efPtr[3] = deta;
  efPtr[4] = phislope;
  efPtr[5] = rphislope;
}

/// Compact the candidate edges selected by @p mask into the final
/// [2, nEdges] int64 edge index and compute their edge features in the same
/// pass. The compaction preserves the order of the candidates.
template <typename T>
__global__ void compactEdgesAndMakeFeatures(
    std::size_t nCandidates, const bool *__restrict__ mask,
    const int *__restrict__ maskSum, const int *__restrict__ srcCandidates,
    const int *__restrict__ tgtCandidates, std::size_t nEdges,
    std::size_t nNodeFeatures, const T *__restrict__ nodeFeatures,
    std::int64_t *edgeIndex, T *edgeFeatures) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nCandidates || !mask[i]) {
    return;
  }

  const std::size_t out = maskSum[i];
  const int src = srcCandidates[i];
  const int tgt = tgtCandidates[i];

  edgeIndex[out] = src;
  edgeIndex[nEdges + out] = tgt;
  computeEdgeFeatures(src, tgt, nNodeFeatures, nodeFeatures,
                      edgeFeatures + out * g_nEdgeFeatures);
}

template <typename T>
__global__ void preprocessHitFeatures(std::size_t nbHits,
                                      std::size_t nNodeFeatures,
                                      const T *nodeFeatures, T *cuda_R,
                                      T *cuda_phi, T *cuda_z, T *cuda_eta,
                                      T *cuda_x, T *cuda_y, T rScale,
                                      T phiScale, T zScale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  enum NodeFeatures { r = 0, phi, z, eta };
  const T *hitFeatures = nodeFeatures + i * nNodeFeatures;

  const T scaledR = hitFeatures[r] * rScale;
  const T scaledPhi = hitFeatures[phi] * phiScale;
  const T scaledZ = hitFeatures[z] * zScale;

  cuda_R[i] = scaledR;
  cuda_phi[i] = scaledPhi;
  cuda_z[i] = scaledZ;
  cuda_eta[i] = hitFeatures[eta];

  const double rd = scaledR;
  const double phid = scaledPhi;
  cuda_x[i] = static_cast<T>(rd * std::cos(phid));
  cuda_y[i] = static_cast<T>(rd * std::sin(phid));
}

static __global__ void mapModuleIdsToNbHits(int *nbHitsOnModule,
                                            std::size_t nHits,
                                            const std::uint64_t *moduleIds,
                                            std::size_t moduleMapSize,
                                            const std::uint64_t *moduleMapKey,
                                            const int *moduleMapVal) {
  auto i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nHits) {
    return;
  }

  auto mId = moduleIds[i];

  int left = 0;
  int right = moduleMapSize - 1;
  while (left <= right) {
    int mid = left + (right - left) / 2;
    if (moduleMapKey[mid] == mId) {
      atomicAdd(&nbHitsOnModule[moduleMapVal[mid]], 1);
      return;
    }
    if (moduleMapKey[mid] < mId) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
}

}  // namespace ActsPlugins::detail
