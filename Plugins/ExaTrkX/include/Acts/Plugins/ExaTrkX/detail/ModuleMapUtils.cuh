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

#define NEW_OPTIMIZATIONS

namespace Acts::detail {

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

inline __global__ void setHitId(std::size_t nHits, std::uint64_t *hit_ids) {
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
  cudaMemcpy(data_cpu.data(), data, size * sizeof(T), cudaMemcpyDeviceToHost);
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
#ifdef NEW_OPTIMIZATIONS
__launch_bounds__(512, 2)
#endif
    count_doublet_edges(int nb_doublets, const int *modules1,
                        const int *modules2, const T *R, const T *z,
                        const T *eta, const T *phi, T *z0_min, T *z0_max,
                        T *deta_min, T *deta_max, T *phi_slope_min,
                        T *phi_slope_max, T *dphi_min, T *dphi_max,
                        const int *indices, T pi, int *nb_edges_total,
                        int *nb_edges_doublet) {
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
                             phi[l], pi, z0, phi_slope, deta, dphi);

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
#ifdef NEW_OPTIMIZATIONS
__launch_bounds__(512, 2)
#endif
    doublet_cuts_new(int nb_doublets, const int *modules1, const int *modules2,
                     const T *R, const T *z, const T *eta, const T *phi,
                     T *z0_min, T *z0_max, T *deta_min, T *deta_max,
                     T *phi_slope_min, T *phi_slope_max, T *dphi_min,
                     T *dphi_max, const int *indices, T pi, int *M1_SP,
                     int *M2_SP, const int *nb_edges) {
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
                             phi[l], pi, z0, phi_slope, deta, dphi);

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
__global__ void
#ifdef NEW_OPTIMIZATIONS
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
                     int *edge_indices, bool *vertices, bool *edge_tag)
//---------------------------------------------------------------------------------------------------------------------------------------------------
{
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
    int p = sorted_M2_SP[k];
    int SP1 = M1_SP[p];
    int SP2 = M2_SP[p];
    bool next_ind = false;
    if (k < last12) {
      next_ind = (SP2 != (M2_SP[sorted_M2_SP[k + 1]]));
    }

    if (!apply_geometric_cuts(
            i, z0[p], phi_slope[p], deta[p], dphi[p], MD12_z0_min, MD12_z0_max,
            MD12_deta_min, MD12_deta_max, MD12_phi_slope_min,
            MD12_phi_slope_max, MD12_dphi_min, MD12_dphi_max)) {
      continue;
    }

    int l = shift23;
#ifdef NEW_OPTIMIZATIONS
    findFirstWithBisect(shift23, ind23 - 1, SP2, l, M1_SP);
#else
    for (; l < ind23 && SP2 != M1_SP[l]; l++) {
    }  // search first hit indice on
#endif

    bool new_elt = false;
    for (; l < ind23 && SP2 == M1_SP[l]; l++) {
      int SP3 = M2_SP[l];
      if (!apply_geometric_cuts(
              i, z0[l], phi_slope[l], deta[l], dphi[l], MD23_z0_min,
              MD23_z0_max, MD23_deta_min, MD23_deta_max, MD23_phi_slope_min,
              MD23_phi_slope_max, MD23_dphi_min, MD23_dphi_max)) {
        continue;
      }

      T diff_dydx = Diff_dydx(x, y, z, SP1, SP2, SP3);
      if (!((diff_dydx >= diff_dydx_min[i]) *
            (diff_dydx <= diff_dydx_max[i]))) {
        continue;
      }

      T diff_dzdr = Diff_dzdr(R, z, SP1, SP2, SP3);
      if (!((diff_dzdr >= diff_dzdr_min[i]) *
            (diff_dzdr <= diff_dzdr_max[i]))) {
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

template <typename T, typename F>
__device__ void doublet_cut_kernel(int i, int sum_nb_src_hits_per_doublet,
                                   int nb_doublets, const int *doublet_offsets,
                                   const int *modules1, const int *modules2,
                                   const T *R, const T *z, const T *eta,
                                   const T *phi, T *z0_min, T *z0_max,
                                   T *deta_min, T *deta_max, T *phi_slope_min,
                                   T *phi_slope_max, T *dphi_min, T *dphi_max,
                                   const int *indices, F &&function) {
  // Since doublet_offsets should be the prefix sum of the
  // number of space points per doublet, we can construct the doublet index
  // TODO this is a linear search, can be optimized
  int doublet_idx = 0;
#if 0
  for (; ; doublet_idx++) {
    if( i < doublet_offsets[doublet_idx+1] ) {
      break;
    }
  }
#else
  {
    // do search with bisection
    int left = 0;
    int right = nb_doublets;
    while (left <= right) {
      int mid = left + (right - left) / 2;
      if (i >= doublet_offsets[mid] && i < doublet_offsets[mid + 1]) {
        doublet_idx = mid;
        break;
      }
      if (i < doublet_offsets[mid]) {
        right = mid - 1;
      } else {
        left = mid + 1;
      }
    }
  }
#endif

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
                           phi[l], detail::g_pi, z0, phi_slope, deta, dphi);

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
    T *dphi_max, const int *indices, int *nb_edges_per_src_hit) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sum_nb_src_hits_per_doublet) {
    return;
  }

  int edges = 0;
  doublet_cut_kernel<T>(i, sum_nb_src_hits_per_doublet, nb_doublets,
                        doublet_offsets, modules1, modules2, R, z, eta, phi,
                        z0_min, z0_max, deta_min, deta_max, phi_slope_min,
                        phi_slope_max, dphi_min, dphi_max, indices,
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
    int *edge_sum) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sum_nb_src_hits_per_doublet) {
    return;
  }

  int edges = 0;
  doublet_cut_kernel<T>(i, sum_nb_src_hits_per_doublet, nb_doublets,
                        doublet_offsets, modules1, modules2, R, z, eta, phi,
                        z0_min, z0_max, deta_min, deta_max, phi_slope_min,
                        phi_slope_max, dphi_min, dphi_max, indices,
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

}  // namespace Acts::detail
