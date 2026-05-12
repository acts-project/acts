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
#include <numbers>
#include <vector>

#include <MMG/CUDA_graph_creator>
#include <cuda_runtime_api.h>

#define USE_LAUNCH_BOUNDS

namespace ActsPlugins::detail {

constexpr float g_pi = std::numbers::pi_v<float>;

template <class T>
__global__ void rescaleFeature(std::size_t nbHits, T *data, T scale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  data[i] *= scale;
}

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

inline void __global__ mapModuleIdsToNbHits(int *nbHitsOnModule,
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

  // bisect moduleMapKey to find mId
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

// ================================
// Triplet cuts (uses sorted_M2_SP)
// ================================

template <typename T>
__device__ void triplet_cuts_inner_loop_body(
    int i, int k, int last12, int shift23, int ind23, T *x, T *y, T *z, T *R,
    T *z0, T *phi_slope, T *deta, T *dphi, T *MD12_z0_min, T *MD12_z0_max,
    T *MD12_deta_min, T *MD12_deta_max, T *MD12_phi_slope_min,
    T *MD12_phi_slope_max, T *MD12_dphi_min, T *MD12_dphi_max, T *MD23_z0_min,
    T *MD23_z0_max, T *MD23_deta_min, T *MD23_deta_max, T *MD23_phi_slope_min,
    T *MD23_phi_slope_max, T *MD23_dphi_min, T *MD23_dphi_max, T *diff_dydx_min,
    T *diff_dydx_max, T *diff_dzdr_min, T *diff_dzdr_max, T pi, int *M1_SP,
    int *M2_SP, int *sorted_M2_SP, int *edge_indices, bool *vertices,
    bool *edge_tag, T epsilon) {
  int p = sorted_M2_SP[k];

  int SP1 = M1_SP[p];
  int SP2 = M2_SP[p];
  bool next_ind = false;
  if (k < last12) {
    next_ind = (SP2 != (M2_SP[sorted_M2_SP[k + 1]]));
  }

  if (!apply_geometric_cuts(i, z0[p], phi_slope[p], deta[p], dphi[p],
                            MD12_z0_min, MD12_phi_slope_min, MD12_deta_min,
                            MD12_dphi_min, MD12_z0_max, MD12_phi_slope_max,
                            MD12_deta_max, MD12_dphi_max)) {
    return;
  }

  int l = shift23;
  findFirstWithBisect(shift23, ind23 - 1, SP2, l, M1_SP);

  bool new_elt = false;
  for (; l < ind23 && SP2 == M1_SP[l]; l++) {
    int SP3 = M2_SP[l];
    if (!apply_geometric_cuts(i, z0[l], phi_slope[l], deta[l], dphi[l],
                              MD23_z0_min, MD23_phi_slope_min, MD23_deta_min,
                              MD23_dphi_min, MD23_z0_max, MD23_phi_slope_max,
                              MD23_deta_max, MD23_dphi_max)) {
      continue;
    }

    T diff_dydx = Diff_dydx(x, y, z, SP1, SP2, SP3, epsilon);
    if (!((diff_dydx >= diff_dydx_min[i]) * (diff_dydx <= diff_dydx_max[i]))) {
      continue;
    }

    T diff_dzdr = Diff_dzdr(R, z, SP1, SP2, SP3, epsilon);
    if (!((diff_dzdr >= diff_dzdr_min[i]) * (diff_dzdr <= diff_dzdr_max[i]))) {
      continue;
    }

    vertices[SP3] = edge_tag[l] = true;
    new_elt = true;
  }
  if (new_elt) {
    edge_tag[p] = vertices[SP1] = vertices[SP2] = true;
  }
  if (next_ind && new_elt) {
    shift23 = l;
  }
}

template <typename T>
__global__ void
#ifdef USE_LAUNCH_BOUNDS
__launch_bounds__(512, 2)
#endif
    triplet_cuts(int nb_src_hits_per_triplet_sum, int nb_triplets,
                 const int *triplet_offsets, const int *modules12_map,
                 const int *modules23_map, T *x, T *y, T *z, T *R, T *z0,
                 T *phi_slope, T *deta, T *dphi, T *MD12_z0_min, T *MD12_z0_max,
                 T *MD12_deta_min, T *MD12_deta_max, T *MD12_phi_slope_min,
                 T *MD12_phi_slope_max, T *MD12_dphi_min, T *MD12_dphi_max,
                 T *MD23_z0_min, T *MD23_z0_max, T *MD23_deta_min,
                 T *MD23_deta_max, T *MD23_phi_slope_min, T *MD23_phi_slope_max,
                 T *MD23_dphi_min, T *MD23_dphi_max, T *diff_dydx_min,
                 T *diff_dydx_max, T *diff_dzdr_min, T *diff_dzdr_max, T pi,
                 int *M1_SP, int *M2_SP, int *sorted_M2_SP, int *edge_indices,
                 bool *vertices, bool *edge_tag, T epsilon) {
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (ii >= nb_src_hits_per_triplet_sum) {
    return;
  }

  // Find triplet index
  int triplet_index = 0;
  locateInPrefixSumBisect(0, nb_triplets, ii, triplet_index, triplet_offsets);

  int module12 = modules12_map[triplet_index];
  int module23 = modules23_map[triplet_index];
  int nb_hits_M12 = edge_indices[module12 + 1] - edge_indices[module12];
  int nb_hits_M23 = edge_indices[module23 + 1] - edge_indices[module23];

  bool hits_on_modules = nb_hits_M12 * nb_hits_M23;
  if (!hits_on_modules) {
    return;
  }

  int shift12 = edge_indices[module12];
  int shift23 = edge_indices[module23];

  int last12 = shift12 + nb_hits_M12 - 1;
  int ind23 = shift23 + nb_hits_M23;

  const int k = shift12 + (ii - triplet_offsets[triplet_index]);

  triplet_cuts_inner_loop_body(
      triplet_index, k, last12, shift23, ind23, x, y, z, R, z0, phi_slope, deta,
      dphi, MD12_z0_min, MD12_z0_max, MD12_deta_min, MD12_deta_max,
      MD12_phi_slope_min, MD12_phi_slope_max, MD12_dphi_min, MD12_dphi_max,
      MD23_z0_min, MD23_z0_max, MD23_deta_min, MD23_deta_max,
      MD23_phi_slope_min, MD23_phi_slope_max, MD23_dphi_min, MD23_dphi_max,
      diff_dydx_min, diff_dydx_max, diff_dzdr_min, diff_dzdr_max, pi, M1_SP,
      M2_SP, sorted_M2_SP, edge_indices, vertices, edge_tag, epsilon);
}

}  // namespace ActsPlugins::detail
