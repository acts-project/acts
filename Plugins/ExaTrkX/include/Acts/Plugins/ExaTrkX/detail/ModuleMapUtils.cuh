// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <CUDA_graph_creator>
#include <algorithm>
#include <iostream>
#include <vector>

#include <cuda_runtime_api.h>

#define USE_LAUNCH_BOUNDS

namespace Acts::detail {

template <class T>
__global__ void rescaleFeature(std::size_t nbHits, T *data, T scale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  data[i] *= scale;
}

constexpr float g_pi = std::numbers::pi_v<float>;

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
      // atomic add 1 to hitIndice[moduleMapVal[mid]]
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

/// Counting kernel to allow counting the edges
template <class T>
__global__ void
#ifdef USE_LAUNCH_BOUNDS
__launch_bounds__(512, 2)
#endif
    count_doublet_edges(int nb_doublets, const int *modules1,
                        const int *modules2, const T *R, const T *z,
                        const T *eta, const T *phi, T *z0_min, T *z0_max,
                        T *deta_min, T *deta_max, T *phi_slope_min,
                        T *phi_slope_max, T *dphi_min, T *dphi_max,
                        const int *indices, T pi, int *nb_edges_total,
                        int *nb_edges_doublet, T epsilon) {
  // loop over module1 SP
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_doublets) {
    return;
  }

  int module1 = modules1[i];
  int module2 = modules2[i];
  int edges = 0;

  for (int k = indices[module1]; k < indices[module1 + 1]; k++) {
    T phi_SP1 = phi[k];
    T eta_SP1 = eta[k];
    T R_SP1 = R[k];
    T z_SP1 = z[k];

    for (int l = indices[module2]; l < indices[module2 + 1]; l++) {
      T z0, phi_slope, deta, dphi;
      hits_geometric_cuts<T>(R_SP1, R[l], z_SP1, z[l], eta_SP1, eta[l], phi_SP1,
                             phi[l], pi, z0, phi_slope, deta, dphi, epsilon);

      if (apply_geometric_cuts(i, z0, phi_slope, deta, dphi, z0_min, z0_max,
                               deta_min, deta_max, phi_slope_min, phi_slope_max,
                               dphi_min, dphi_max)) {
        edges++;
      }
    }
  }

  // increase global and local counter
  nb_edges_doublet[i] = edges;
  atomicAdd(nb_edges_total, edges);
}

/// New kernel that use precounted number of edges
template <class T>
__global__ void
#ifdef USE_LAUNCH_BOUNDS
__launch_bounds__(512, 2)
#endif
    doublet_cuts_new(int nb_doublets, const int *modules1, const int *modules2,
                     const T *R, const T *z, const T *eta, const T *phi,
                     T *z0_min, T *z0_max, T *deta_min, T *deta_max,
                     T *phi_slope_min, T *phi_slope_max, T *dphi_min,
                     T *dphi_max, const int *indices, T pi, int *M1_SP,
                     int *M2_SP, const int *nb_edges, T epsilon) {
  // loop over module1 SP
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_doublets) {
    return;
  }

  int module1 = modules1[i];
  int module2 = modules2[i];
  int edges = nb_edges[i];

  for (int k = indices[module1]; k < indices[module1 + 1]; k++) {
    T phi_SP1 = phi[k];
    T eta_SP1 = eta[k];
    T R_SP1 = R[k];
    T z_SP1 = z[k];

    for (int l = indices[module2]; l < indices[module2 + 1]; l++) {
      T z0, phi_slope, deta, dphi;
      hits_geometric_cuts<T>(R_SP1, R[l], z_SP1, z[l], eta_SP1, eta[l], phi_SP1,
                             phi[l], pi, z0, phi_slope, deta, dphi, epsilon);

      if (apply_geometric_cuts(i, z0, phi_slope, deta, dphi, z0_min, z0_max,
                               deta_min, deta_max, phi_slope_min, phi_slope_max,
                               dphi_min, dphi_max)) {
        M1_SP[edges] = k;
        M2_SP[edges] = l;
        edges++;
      }
    }
  }
}

__device__ void findFirstWithBisect(int left, int right, int query, int &result,
                                    const int *array) {
  while (left <= right) {
    int mid = left + (right - left) / 2;
    // only terminate search if we found the first index of the hit
    // guard against array[mid-1] to be outside of the array
    if (array[mid] == query && (mid == left || array[mid - 1] != query)) {
      result = mid;
      break;
    }
    if (array[mid] < query) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
}

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
                            MD12_z0_min, MD12_z0_max, MD12_deta_min,
                            MD12_deta_max, MD12_phi_slope_min,
                            MD12_phi_slope_max, MD12_dphi_min, MD12_dphi_max)) {
    return;
  }

  int l = shift23;
#ifdef USE_LAUNCH_BOUNDS
  findFirstWithBisect(shift23, ind23 - 1, SP2, l, M1_SP);
#else
  for (; l < ind23 && SP2 != M1_SP[l]; l++) {
  }  // search first hit indice on
#endif

  bool new_elt = false;
  for (; l < ind23 && SP2 == M1_SP[l]; l++) {
    int SP3 = M2_SP[l];
    if (!apply_geometric_cuts(
            i, z0[l], phi_slope[l], deta[l], dphi[l], MD23_z0_min, MD23_z0_max,
            MD23_deta_min, MD23_deta_max, MD23_phi_slope_min,
            MD23_phi_slope_max, MD23_dphi_min, MD23_dphi_max)) {
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
    triplet_cuts_new(int nb_triplets, const int *modules12_map,
                     const int *modules23_map, T *x, T *y, T *z, T *R, T *z0,
                     T *phi_slope, T *deta, T *dphi, T *MD12_z0_min,
                     T *MD12_z0_max, T *MD12_deta_min, T *MD12_deta_max,
                     T *MD12_phi_slope_min, T *MD12_phi_slope_max,
                     T *MD12_dphi_min, T *MD12_dphi_max, T *MD23_z0_min,
                     T *MD23_z0_max, T *MD23_deta_min, T *MD23_deta_max,
                     T *MD23_phi_slope_min, T *MD23_phi_slope_max,
                     T *MD23_dphi_min, T *MD23_dphi_max, T *diff_dydx_min,
                     T *diff_dydx_max, T *diff_dzdr_min, T *diff_dzdr_max, T pi,
                     int *M1_SP, int *M2_SP, int *sorted_M2_SP,
                     int *edge_indices, bool *vertices, bool *edge_tag,
                     T epsilon) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_triplets) {
    return;
  }

  int module12 = modules12_map[i];
  int module23 = modules23_map[i];

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

  for (int k = shift12; k <= last12; k++) {
    triplet_cuts_inner_loop_body(
        i, k, last12, shift23, ind23, x, y, z, R, z0, phi_slope, deta, dphi,
        MD12_z0_min, MD12_z0_max, MD12_deta_min, MD12_deta_max,
        MD12_phi_slope_min, MD12_phi_slope_max, MD12_dphi_min, MD12_dphi_max,
        MD23_z0_min, MD23_z0_max, MD23_deta_min, MD23_deta_max,
        MD23_phi_slope_min, MD23_phi_slope_max, MD23_dphi_min, MD23_dphi_max,
        diff_dydx_min, diff_dydx_max, diff_dzdr_min, diff_dzdr_max, pi, M1_SP,
        M2_SP, sorted_M2_SP, edge_indices, vertices, edge_tag, epsilon);
  }
}

// ============================
// Utiles for new doublet edges
// ============================

__global__ void count_src_hits_per_doublet(int nb_doublets, const int *modules1,
                                           const int *indices,
                                           int *nb_src_hits_per_doublet) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_doublets) {
    return;
  }

  int module1 = modules1[i];
  nb_src_hits_per_doublet[i] = indices[module1 + 1] - indices[module1];
}

/// Assume we have an query value i in the range (left, right),
/// and a prefix sum, also in this range.
/// This function finds the lower bound idx of the step in the prefix sum that
/// contains i
__device__ void locateInPrefixSumBisect(int left, int right, int i, int &result,
                                        const int *prefix_sum) {
  while (left <= right) {
    int mid = left + (right - left) / 2;
    if (i >= prefix_sum[mid] && i < prefix_sum[mid + 1]) {
      result = mid;
      break;
    }
    if (i < prefix_sum[mid]) {
      right = mid - 1;
    } else {
      left = mid + 1;
    }
  }
}

template <typename T, typename F>
__device__ void doublet_cut_kernel(
    int i, int sum_nb_src_hits_per_doublet, int nb_doublets,
    const int *doublet_offsets, const int *modules1, const int *modules2,
    const T *R, const T *z, const T *eta, const T *phi, T *z0_min, T *z0_max,
    T *deta_min, T *deta_max, T *phi_slope_min, T *phi_slope_max, T *dphi_min,
    T *dphi_max, const int *indices, T epsilon, F &&function) {
  // Since doublet_offsets should be the prefix sum of the
  // number of space points per doublet, we can construct the doublet index
  // TODO this is a linear search, can be optimized
  int doublet_idx = 0;
  locateInPrefixSumBisect(0, nb_doublets, i, doublet_idx, doublet_offsets);

  const int module1 = modules1[doublet_idx];
  const int module2 = modules2[doublet_idx];

  // we can reconstruct the SP1 index from start-index of SPs for module1,
  // and the difference between the current i and the start of the doublet
  // Note that several doublets can have the same module1
  const int k = indices[module1] + (i - doublet_offsets[doublet_idx]);

  T phi_SP1 = phi[k];
  T eta_SP1 = eta[k];
  T R_SP1 = R[k];
  T z_SP1 = z[k];

  for (int l = indices[module2]; l < indices[module2 + 1]; l++) {
    T z0, phi_slope, deta, dphi;
    hits_geometric_cuts<T>(R_SP1, R[l], z_SP1, z[l], eta_SP1, eta[l], phi_SP1,
                           phi[l], detail::g_pi, z0, phi_slope, deta, dphi,
                           epsilon);

    if (apply_geometric_cuts(doublet_idx, z0, phi_slope, deta, dphi, z0_min,
                             z0_max, deta_min, deta_max, phi_slope_min,
                             phi_slope_max, dphi_min, dphi_max)) {
      function(k, l);
    }
  }
}

template <class T>
__global__ void count_doublet_edges_new(
    int sum_nb_src_hits_per_doublet, int nb_doublets,
    const int *doublet_offsets, const int *modules1, const int *modules2,
    const T *R, const T *z, const T *eta, const T *phi, T *z0_min, T *z0_max,
    T *deta_min, T *deta_max, T *phi_slope_min, T *phi_slope_max, T *dphi_min,
    T *dphi_max, const int *indices, int *nb_edges_per_src_hit, T epsilon) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sum_nb_src_hits_per_doublet) {
    return;
  }

  int edges = 0;
  doublet_cut_kernel<T>(i, sum_nb_src_hits_per_doublet, nb_doublets,
                        doublet_offsets, modules1, modules2, R, z, eta, phi,
                        z0_min, z0_max, deta_min, deta_max, phi_slope_min,
                        phi_slope_max, dphi_min, dphi_max, indices, epsilon,
                        [&] __device__(int, int) { edges++; });
  nb_edges_per_src_hit[i] = edges;
}

template <class T>
__global__ void build_doublet_edges_new(
    int sum_nb_src_hits_per_doublet, int nb_doublets,
    const int *doublet_offsets, const int *modules1, const int *modules2,
    const T *R, const T *z, const T *eta, const T *phi, T *z0_min, T *z0_max,
    T *deta_min, T *deta_max, T *phi_slope_min, T *phi_slope_max, T *dphi_min,
    T *dphi_max, const int *indices, int *reduced_M1_hits, int *reduced_M2_hits,
    int *edge_sum, T epsilon) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sum_nb_src_hits_per_doublet) {
    return;
  }

  int edges = 0;
  doublet_cut_kernel<T>(i, sum_nb_src_hits_per_doublet, nb_doublets,
                        doublet_offsets, modules1, modules2, R, z, eta, phi,
                        z0_min, z0_max, deta_min, deta_max, phi_slope_min,
                        phi_slope_max, dphi_min, dphi_max, indices, epsilon,
                        [&] __device__(int k, int l) {
                          reduced_M1_hits[edge_sum[i] + edges] = k;
                          reduced_M2_hits[edge_sum[i] + edges] = l;
                          edges++;
                        });
}

__global__ void computeDoubletEdgeSum(int nb_doublets,
                                      const int *doublet_offsets,
                                      const int *nb_edges_per_src_hit,
                                      int *edge_sum) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  // Note that we want to get the maximum element es well
  if (i >= nb_doublets + 1) {
    return;
  }

  // Since nb_edges_per_src_hit is already a prefix sum, we can just copy the
  // elements on the boundary positions
  edge_sum[i] = nb_edges_per_src_hit[doublet_offsets[i]];
}

// ================================
// New kernels for the triplet cuts
// ================================

__global__ void count_triplet_hits(int nb_triplets, const int *modules12_map,
                                   const int *modules23_map,
                                   const int *edge_indices,
                                   int *src_hits_per_triplet) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_triplets) {
    return;
  }
  int module12 = modules12_map[i];
  int module23 = modules23_map[i];

  int nb_hits_M12 = edge_indices[module12 + 1] - edge_indices[module12];
  int nb_hits_M23 = edge_indices[module23 + 1] - edge_indices[module23];

  if (nb_hits_M12 == 0 || nb_hits_M23 == 0) {
    src_hits_per_triplet[i] = 0;
    return;
  }

  int shift12 = edge_indices[module12];
  int last12 = shift12 + nb_hits_M12 - 1;

  src_hits_per_triplet[i] = last12 - shift12 + 1;
}

template <typename T>
__global__ void
#ifdef USE_LAUNCH_BOUNDS
__launch_bounds__(512, 2)
#endif
    triplet_cuts_new2(int nb_src_hits_per_triplet_sum, int nb_triplets,
                      const int *triplet_offsets, const int *modules12_map,
                      const int *modules23_map, T *x, T *y, T *z, T *R, T *z0,
                      T *phi_slope, T *deta, T *dphi, T *MD12_z0_min,
                      T *MD12_z0_max, T *MD12_deta_min, T *MD12_deta_max,
                      T *MD12_phi_slope_min, T *MD12_phi_slope_max,
                      T *MD12_dphi_min, T *MD12_dphi_max, T *MD23_z0_min,
                      T *MD23_z0_max, T *MD23_deta_min, T *MD23_deta_max,
                      T *MD23_phi_slope_min, T *MD23_phi_slope_max,
                      T *MD23_dphi_min, T *MD23_dphi_max, T *diff_dydx_min,
                      T *diff_dydx_max, T *diff_dzdr_min, T *diff_dzdr_max,
                      T pi, int *M1_SP, int *M2_SP, int *sorted_M2_SP,
                      int *edge_indices, bool *vertices, bool *edge_tag,
                      T epsilon) {
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

  // TODO does this work???
  const int k = shift12 + (ii - triplet_offsets[triplet_index]);

  // From here on starts the original kernel, with the loop over k pulled out of
  // the kernel
  triplet_cuts_inner_loop_body(
      triplet_index, k, last12, shift23, ind23, x, y, z, R, z0, phi_slope, deta,
      dphi, MD12_z0_min, MD12_z0_max, MD12_deta_min, MD12_deta_max,
      MD12_phi_slope_min, MD12_phi_slope_max, MD12_dphi_min, MD12_dphi_max,
      MD23_z0_min, MD23_z0_max, MD23_deta_min, MD23_deta_max,
      MD23_phi_slope_min, MD23_phi_slope_max, MD23_dphi_min, MD23_dphi_max,
      diff_dydx_min, diff_dydx_max, diff_dzdr_min, diff_dzdr_max, pi, M1_SP,
      M2_SP, sorted_M2_SP, edge_indices, vertices, edge_tag, epsilon);
}

// =======================
// New kernels for sorting
// =======================

__global__ void __launch_bounds__(512, 4)
    block_odd_even_sort(const int *M2_hits, const int *cuda_edge_sum,
                        int *M2_idxs) {
  bool sorted{};
  auto comparison = thrust::less<int>{};

  const int begin = cuda_edge_sum[blockIdx.x];
  const int end = cuda_edge_sum[blockIdx.x + 1];
  if (end - begin == 0) {
    return;
  }

  do {
    sorted = true;

    for (std::uint32_t j =
             begin + 2 * static_cast<std::uint32_t>(threadIdx.x) + 1;
         j < end - 1; j += 2 * blockDim.x) {
      if (comparison(M2_hits[M2_idxs[j + 1]], M2_hits[M2_idxs[j]])) {
        swap(M2_idxs[j + 1], M2_idxs[j]);
        sorted = false;
      }
    }

    __syncthreads();

    for (std::uint32_t j = begin + 2 * static_cast<std::uint32_t>(threadIdx.x);
         j < end - 1; j += 2 * blockDim.x) {
      if (comparison(M2_hits[M2_idxs[j + 1]], M2_hits[M2_idxs[j]])) {
        swap(M2_idxs[j + 1], M2_idxs[j]);
        sorted = false;
      }
    }
  } while (__syncthreads_or(!sorted));
}

}  // namespace Acts::detail
