// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

#include <cuda_runtime_api.h>

/// CUDA kernels for the module map graph construction in ModuleMapCuda.
namespace ActsPlugins::detail {

// Precision strategy
// ------------------
// The cut values come from the reference double-precision implementation, but
// consumer GPUs run FP64 at a small fraction of the FP32 rate. The kernels
// below therefore evaluate everything in float, in a way that reproduces the
// double-precision decisions exactly:
//  - Quantities whose float evaluation is exact, or equals the correctly
//    rounded double result (single float divisions of float inputs), are
//    compared directly. This requires IEEE division (no --use_fast_math).
//  - Other quantities are float approximations with a proven error bound.
//    A cut is decided from the approximation unless it lies within that bound
//    of a cut edge; only those rare cases are evaluated in double precision.

// Maximum number of threads per block of the graph construction kernels
// (their launch bound), i.e. the maximum of ModuleMapCuda::Config::gpuBlocks
constexpr int kMaxBlockSize = 512;

// Relative error margin for the float approximations (12 float ulps). The
// error bounds it has to cover are given next to each use.
constexpr float kFilterErrScale = 6.f * std::numeric_limits<float>::epsilon();

// Reference double-precision geometry of one doublet edge, used when a float
// approximation is too close to a cut edge. Must stay arithmetically identical
// to the reference implementation.
template <class T>
struct EdgeGeoReference {
  T z0;
  T phiSlope;
  double dydx;
  double dzdr;
};

template <class T>
__device__ EdgeGeoReference<T> edgeGeoReference(
    int SP1, int SP2, const T *__restrict__ R, const T *__restrict__ z,
    const T *__restrict__ x, const T *__restrict__ y, const T *__restrict__ phi,
    double pi, T epsilon) {
  double R1 = R[SP1];
  T R2 = R[SP2];
  T z1 = z[SP1];
  T z2 = z[SP2];

  double dr = R2 - R1;
  double dz = z2 - z1;
  double dphiLoc = phi[SP2] - phi[SP1];
  T pipi = 2 * pi;

  if (dphiLoc > pi) {
    dphiLoc -= pipi;
  } else if (dphiLoc < -pi) {
    dphiLoc += pipi;
  }

  EdgeGeoReference<T> ref{};
  if (abs(dr) > epsilon) {
    ref.phiSlope = dphiLoc / dr;
    ref.z0 = z1 - R1 * dz / dr;
  } else {
    ref.z0 = 0;
    ref.phiSlope = 0;
  }

  double dx = x[SP2] - x[SP1];
  const double dy = y[SP2] - y[SP1];
  double drSlope = R[SP2] - R[SP1];
  const double dzSlope = z[SP2] - z[SP1];
  if (abs(dx) < epsilon) {
    dx = 0;
  }
  if (abs(drSlope) < epsilon) {
    drSlope = 0;
  }
  ref.dydx = dx ? dy / dx : 0;
  ref.dzdr = drSlope ? dzSlope / drSlope : 0;
  return ref;
}

// Float geometry of each doublet edge, used by the triplet cuts:
//  - geo:   (z0, phiSlope, deta, dphi); deta and dphi are exact, z0 and
//           phiSlope are float approximations of the reference.
//  - slope: (dydx, dzdr, z0 error margin, unused); dydx and dzdr are the
//           correctly rounded float values of the reference slopes.
template <class T>
__global__ __launch_bounds__(kMaxBlockSize, 3) void hitsGeometricCutsPacked(
    float4 *geo, float4 *slope, const int *__restrict__ SPi,
    const int *__restrict__ SPo, const T *__restrict__ R,
    const T *__restrict__ z, const T *__restrict__ x, const T *__restrict__ y,
    const T *__restrict__ eta, const T *__restrict__ phi, T pi, T epsilon,
    int nbDoublets) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nbDoublets) {
    return;
  }

  const int SP1 = SPi[i];
  const int SP2 = SPo[i];
  const T R1 = R[SP1];
  const T R2 = R[SP2];
  const T z1 = z[SP1];
  const T z2 = z[SP2];

  const T deta = eta[SP2] - eta[SP1];
  T dphi = phi[SP2] - phi[SP1];
  if (dphi > pi) {
    dphi -= 2 * pi;
  } else if (dphi < -pi) {
    dphi += 2 * pi;
  }

  const T dr = R2 - R1;
  const T dz = z2 - z1;

  T z0 = 0;
  T phiSlope = 0;
  T z0Err = 0;
  // The reference tests the double difference R2 - R1 against epsilon
  if (abs(static_cast<double>(R2) - static_cast<double>(R1)) > epsilon) {
    // With u = 2^-24: phiSlope is within ~3u|phiSlope| of the reference,
    // z0 within ~5u|q| + 2u|z1|; both are covered by kFilterErrScale
    phiSlope = dphi / dr;
    const T q = R1 * dz / dr;
    z0 = z1 - q;
    z0Err = kFilterErrScale * (fabs(z1) + fabs(q)) + T{1e-30};
  }
  geo[i] = make_float4(z0, phiSlope, deta, dphi);

  // dx, dy, dr and dz are float differences in the reference as well
  const T dx = x[SP2] - x[SP1];
  const T dy = y[SP2] - y[SP1];
  const T dydx = (fabs(dx) < epsilon) ? T{0} : dy / dx;
  const T dzdr = (fabs(dr) < epsilon) ? T{0} : dz / dr;
  slope[i] = make_float4(dydx, dzdr, z0Err, 0.f);
}

static __global__ void countSourceHitsPerDoublet(int *nbSrcHitsPerDoublet,
                                                 const int *modules1,
                                                 const int *indices,
                                                 int nbDoublets) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nbDoublets) {
    return;
  }

  int module1 = modules1[i];
  nbSrcHitsPerDoublet[i] = indices[module1 + 1] - indices[module1];
}

static __global__ void compactActiveDoublets(int *activeDoublets,
                                             const int *activeFlags,
                                             const int *activeOffsets,
                                             int nbDoublets) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nbDoublets || activeFlags[i] == 0) {
    return;
  }

  activeDoublets[activeOffsets[i]] = i;
}

// Fill workToItem[offsets[i]..offsets[i+1]) with i. A group of
// kThreadsPerItem threads cooperates on each item, so the writes are coalesced
// and long ranges do not serialize on a single thread.
template <int kThreadsPerItem>
__global__ void buildWorkToItem(int *workToItem,
                                const int *__restrict__ offsets, int nbItems) {
  const std::size_t tid =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  const std::size_t item = tid / kThreadsPerItem;
  const int lane = tid % kThreadsPerItem;
  if (item >= static_cast<std::size_t>(nbItems)) {
    return;
  }
  const int begin = offsets[item];
  const int end = offsets[item + 1];
  for (int i = begin + lane; i < end; i += kThreadsPerItem) {
    workToItem[i] = static_cast<int>(item);
  }
}

// Doublet pair selection shared by countDoubletEdges and
// buildDoubletEdgesActive. The phi wrap is exact in float (Sterbenz lemma),
// dr is a float difference in the reference too, and phiSlope is a single
// float division; only z0 needs the filtered evaluation.
template <class T>
__device__ __forceinline__ bool doubletPairPasses(
    T R1, T z1, T eta1, T phi1, const T *__restrict__ R,
    const T *__restrict__ z, const T *__restrict__ eta,
    const T *__restrict__ phi, int l, T z0Min, T z0Max, T detaMin, T detaMax,
    T phiSlopeMin, T phiSlopeMax, T dphiMin, T dphiMax, T pi, T epsilon) {
  const T deta = eta[l] - eta1;
  if (!((detaMin <= deta) && (deta <= detaMax))) {
    return false;
  }

  T dphi = phi[l] - phi1;
  if (dphi > pi) {
    dphi -= 2 * pi;
  } else if (dphi < -pi) {
    dphi += 2 * pi;
  }
  if (!((dphiMin <= dphi) && (dphi <= dphiMax))) {
    return false;
  }

  // R and z are only loaded for pairs that pass the eta/phi windows
  const T dr = R[l] - R1;
  if (!(fabs(dr) > epsilon)) {
    // Degenerate dr: z0 and phiSlope are defined as 0
    return (z0Min <= T{0}) && (T{0} <= z0Max) && (phiSlopeMin <= T{0}) &&
           (T{0} <= phiSlopeMax);
  }

  const T phiSlope = dphi / dr;
  if (!((phiSlopeMin <= phiSlope) && (phiSlope <= phiSlopeMax))) {
    return false;
  }

  // With u = 2^-24, the float z0 differs from the rounded double reference by
  // less than 6u(|z1| + |q|); the margin is 8u(|z1| + |q|). Non-finite values
  // fail the margin test and take the double-precision path.
  const T dz = z[l] - z1;
  const T q = R1 * dz / dr;
  const T z0Fast = z1 - q;
  const T margin =
      T{4} * std::numeric_limits<T>::epsilon() * (fabs(z1) + fabs(q)) +
      T{1e-30};
  if (fabs(z0Fast - z0Min) > margin && fabs(z0Fast - z0Max) > margin) {
    return (z0Min <= z0Fast) && (z0Fast <= z0Max);
  }

  const T z0 = z1 - static_cast<double>(R1) * dz / static_cast<double>(dr);
  return (z0Min <= z0) && (z0 <= z0Max);
}

// The count pass records the accepted hit pairs as a bit mask per source hit,
// so the build pass can write the edges without evaluating the cuts again.
// Target modules with more hits than fit in the mask are re-evaluated instead.
constexpr int kMaxMaskedTargetHits = 64;

template <class T>
__global__ __launch_bounds__(kMaxBlockSize, 3) void countDoubletEdges(
    int *nbEdgesPerSrcHit, std::uint64_t *pairMasks,
    const int *__restrict__ srcWorkToDoublet,
    const int *__restrict__ doubletOffsets, const int *__restrict__ modules1,
    const int *__restrict__ modules2, const T *__restrict__ R,
    const T *__restrict__ z, const T *__restrict__ eta,
    const T *__restrict__ phi, const T *__restrict__ z0Min,
    const T *__restrict__ detaMin, const T *__restrict__ phiSlopeMin,
    const T *__restrict__ dphiMin, const T *__restrict__ z0Max,
    const T *__restrict__ detaMax, const T *__restrict__ phiSlopeMax,
    const T *__restrict__ dphiMax, const int *__restrict__ indices, T pi,
    T epsilon, int sumNbSrcHitsPerDoublet) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sumNbSrcHitsPerDoublet) {
    return;
  }

  const int doubletIdx = srcWorkToDoublet[i];
  const int module1 = modules1[doubletIdx];
  const int module2 = modules2[doubletIdx];
  const int srcBegin = doubletOffsets[doubletIdx];
  const int k = indices[module1] + (i - srcBegin);
  const T phi_SP1 = phi[k];
  const T eta_SP1 = eta[k];
  const T R_SP1 = R[k];
  const T z_SP1 = z[k];
  const T z0MinV = z0Min[doubletIdx];
  const T detaMinV = detaMin[doubletIdx];
  const T phiSlopeMinV = phiSlopeMin[doubletIdx];
  const T dphiMinV = dphiMin[doubletIdx];
  const T z0MaxV = z0Max[doubletIdx];
  const T detaMaxV = detaMax[doubletIdx];
  const T phiSlopeMaxV = phiSlopeMax[doubletIdx];
  const T dphiMaxV = dphiMax[doubletIdx];

  const int begin2 = indices[module2];
  const int end2 = indices[module2 + 1];
  const bool useMask = (end2 - begin2) <= kMaxMaskedTargetHits;

  int edges = 0;
  std::uint64_t mask = 0;
  for (int l = begin2; l < end2; l++) {
    const bool pass =
        doubletPairPasses<T>(R_SP1, z_SP1, eta_SP1, phi_SP1, R, z, eta, phi, l,
                             z0MinV, z0MaxV, detaMinV, detaMaxV, phiSlopeMinV,
                             phiSlopeMaxV, dphiMinV, dphiMaxV, pi, epsilon);
    edges += pass;
    if (useMask && pass) {
      mask |= std::uint64_t{1} << (l - begin2);
    }
  }

  nbEdgesPerSrcHit[i] = edges;
  pairMasks[i] = mask;
}

static __global__ void markActiveSrcWork(int *flags,
                                         const int *__restrict__ edgeSum,
                                         int sumNbSrcHitsPerDoublet) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= sumNbSrcHitsPerDoublet) {
    return;
  }
  flags[i] = (edgeSum[i + 1] > edgeSum[i]) ? 1 : 0;
}

template <class T>
__global__ __launch_bounds__(kMaxBlockSize, 3) void buildDoubletEdgesActive(
    int *reduced_M1_hits, int *reduced_M2_hits, int nbActiveSrc,
    const std::uint64_t *__restrict__ pairMasks,
    const int *__restrict__ activeSrcWork,
    const int *__restrict__ srcWorkToDoublet,
    const int *__restrict__ doubletOffsets, const int *__restrict__ indices,
    const int *__restrict__ edgeSum, const int *__restrict__ modules1,
    const int *__restrict__ modules2, const T *__restrict__ R,
    const T *__restrict__ z, const T *__restrict__ eta,
    const T *__restrict__ phi, const T *__restrict__ z0Min,
    const T *__restrict__ detaMin, const T *__restrict__ phiSlopeMin,
    const T *__restrict__ dphiMin, const T *__restrict__ z0Max,
    const T *__restrict__ detaMax, const T *__restrict__ phiSlopeMax,
    const T *__restrict__ dphiMax, T pi, T epsilon) {
  const int a = blockIdx.x * blockDim.x + threadIdx.x;
  if (a >= nbActiveSrc) {
    return;
  }

  const int srcWorkI = activeSrcWork[a];
  const int doubletIdx = srcWorkToDoublet[srcWorkI];
  const int module1 = modules1[doubletIdx];
  const int module2 = modules2[doubletIdx];
  const int k = indices[module1] + (srcWorkI - doubletOffsets[doubletIdx]);
  const int begin2 = indices[module2];
  const int end2 = indices[module2 + 1];
  int out = edgeSum[srcWorkI];

  if ((end2 - begin2) <= kMaxMaskedTargetHits) {
    // Emit the pairs recorded by the count pass, lowest bit first, which is
    // the same order as the loop below
    std::uint64_t mask = pairMasks[srcWorkI];
    while (mask != 0) {
      const int bit = __ffsll(static_cast<long long>(mask)) - 1;
      reduced_M1_hits[out] = k;
      reduced_M2_hits[out] = begin2 + bit;
      out++;
      mask &= mask - 1;
    }
    return;
  }

  const T phi_SP1 = phi[k];
  const T eta_SP1 = eta[k];
  const T R_SP1 = R[k];
  const T z_SP1 = z[k];
  const T z0MinV = z0Min[doubletIdx];
  const T detaMinV = detaMin[doubletIdx];
  const T phiSlopeMinV = phiSlopeMin[doubletIdx];
  const T dphiMinV = dphiMin[doubletIdx];
  const T z0MaxV = z0Max[doubletIdx];
  const T detaMaxV = detaMax[doubletIdx];
  const T phiSlopeMaxV = phiSlopeMax[doubletIdx];
  const T dphiMaxV = dphiMax[doubletIdx];

  for (int l = begin2; l < end2; l++) {
    if (!doubletPairPasses<T>(R_SP1, z_SP1, eta_SP1, phi_SP1, R, z, eta, phi, l,
                              z0MinV, z0MaxV, detaMinV, detaMaxV, phiSlopeMinV,
                              phiSlopeMaxV, dphiMinV, dphiMaxV, pi, epsilon)) {
      continue;
    }
    reduced_M1_hits[out] = k;
    reduced_M2_hits[out] = l;
    out++;
  }
}

static __global__ void doubletEdgeSum(int *edgeSum, const int *doubletOffsets,
                                      const int *nbEdgesPerSrcHit,
                                      int nbDoublets) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nbDoublets + 1) {
    return;
  }
  edgeSum[i] = nbEdgesPerSrcHit[doubletOffsets[i]];
}

static __global__ void countTripletHits(int *srcHitsPerTriplet,
                                        const int *modules12Map,
                                        const int *modules23Map,
                                        const int *edgeIndices,
                                        int nbTriplets) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nbTriplets) {
    return;
  }

  // One work item per M1->M2 edge, if the M2->M3 doublet has edges at all
  const int module12 = modules12Map[i];
  const int module23 = modules23Map[i];
  const int nb_edges_M12 = edgeIndices[module12 + 1] - edgeIndices[module12];
  const int nb_edges_M23 = edgeIndices[module23 + 1] - edgeIndices[module23];
  srcHitsPerTriplet[i] = (nb_edges_M23 == 0) ? 0 : nb_edges_M12;
}

// Decide lo <= v <= hi for a float approximation v whose reference value is
// within err of it. Returns 1 (pass) or 0 (fail) when v is farther than err
// from both edges, so the reference value gives the same decision, and -1
// when the decision needs the reference value. NaN gives -1.
template <class T>
__device__ __forceinline__ int inRangeFiltered(T v, T err, T lo, T hi) {
  if (fabs(v - lo) > err && fabs(v - hi) > err) {
    return (lo <= v) && (v <= hi);
  }
  return -1;
}

template <class T>
__device__ __forceinline__ bool inRange(T v, T lo, T hi) {
  return (lo <= v) && (v <= hi);
}

// A pair of doublet edges (k: M1->M2, l: M2->M3) of a module triplet whose
// cuts could not be decided from the float geometry
struct TripletFallbackPair {
  int tripletIndex;
  int k;
  int l;
};

// Capacity of the fallback queue (typically ~100 pairs per event are queued).
// On overflow the host reruns the triplet cuts with inline fallback.
constexpr int kTripletFallbackCapacity = 1 << 20;

// Triplet cuts on pairs of doublet edges (k: M1->M2, l: M2->M3), evaluated
// on the float geometry from hitsGeometricCutsPacked. Pairs that cannot be
// decided from it are resolved with edgeGeoReference:
//  - kDeferFallback = true: queued for tripletPairCutsFallback, which keeps
//    double-precision code (and its register pressure) out of this kernel.
//  - kDeferFallback = false: resolved inline; used if the queue overflows.
// The slope differences are within ~3u(|slopeK| + |slopeL|) of the
// reference, u = 2^-24.
template <typename T, bool kDeferFallback>
__global__ void tripletPairCutsFusedM23(
    bool *edgeTag, int nbWorkItems, const int *__restrict__ tripletOffsets,
    const int *__restrict__ workToTriplet, const int *__restrict__ modules12Map,
    const int *__restrict__ modules23Map, const float4 *__restrict__ geo,
    const float4 *__restrict__ edgeSlope, const T *__restrict__ MD12_z0_min,
    const T *__restrict__ MD12_phi_slope_min,
    const T *__restrict__ MD12_deta_min, const T *__restrict__ MD12_dphi_min,
    const T *__restrict__ MD12_z0_max, const T *__restrict__ MD12_phi_slope_max,
    const T *__restrict__ MD12_deta_max, const T *__restrict__ MD12_dphi_max,
    const T *__restrict__ MD23_z0_min, const T *__restrict__ MD23_phi_slope_min,
    const T *__restrict__ MD23_deta_min, const T *__restrict__ MD23_dphi_min,
    const T *__restrict__ MD23_z0_max, const T *__restrict__ MD23_phi_slope_max,
    const T *__restrict__ MD23_deta_max, const T *__restrict__ MD23_dphi_max,
    const T *__restrict__ diffDydxMin, const T *__restrict__ diffDydxMax,
    const T *__restrict__ diffDzdrMin, const T *__restrict__ diffDzdrMax,
    const int *__restrict__ M1_SP, const int *__restrict__ M2_SP,
    const int *__restrict__ edgeIndices,
    const int *__restrict__ doubletSrcOffsets,
    const int *__restrict__ doubletModule1, const int *__restrict__ hitIndices,
    const int *__restrict__ edgeSumPerSrcHit, const T *__restrict__ R,
    const T *__restrict__ z, const T *__restrict__ x, const T *__restrict__ y,
    const T *__restrict__ phi, T pi, T epsilon,
    TripletFallbackPair *fallbackPairs, int *fallbackCount,
    int fallbackCapacity) {
  int workI = blockIdx.x * blockDim.x + threadIdx.x;
  if (workI >= nbWorkItems) {
    return;
  }

  // Work items only exist for triplets where both doublets have edges (see
  // countTripletHits), so no emptiness check is needed here
  const int tripletIndex = workToTriplet[workI];
  const int module12 = modules12Map[tripletIndex];
  const int module23 = modules23Map[tripletIndex];
  const int k = edgeIndices[module12] + (workI - tripletOffsets[tripletIndex]);

  // ---- cuts on edge k (module doublet M1-M2) ----
  const float4 gk = geo[k];
  if (!(inRange<T>(gk.z, MD12_deta_min[tripletIndex],
                   MD12_deta_max[tripletIndex]) &&
        inRange<T>(gk.w, MD12_dphi_min[tripletIndex],
                   MD12_dphi_max[tripletIndex]))) {
    return;
  }

  const float4 sk = edgeSlope[k];
  const T z0Min12 = MD12_z0_min[tripletIndex];
  const T z0Max12 = MD12_z0_max[tripletIndex];
  const T psMin12 = MD12_phi_slope_min[tripletIndex];
  const T psMax12 = MD12_phi_slope_max[tripletIndex];
  const int z0OkK = inRangeFiltered<T>(gk.x, sk.z, z0Min12, z0Max12);
  const int psOkK =
      inRangeFiltered<T>(gk.y, kFilterErrScale * fabs(gk.y), psMin12, psMax12);
  if (z0OkK == 0 || psOkK == 0) {
    return;
  }
  bool refKUndecided = (z0OkK < 0 || psOkK < 0);

  // Inline mode only: the reference geometry of edge k, computed at most once
  EdgeGeoReference<T> refK{};
  bool haveRefK = false;
  if constexpr (!kDeferFallback) {
    if (refKUndecided) {
      refK =
          edgeGeoReference<T>(M1_SP[k], M2_SP[k], R, z, x, y, phi, pi, epsilon);
      haveRefK = true;
      if (!(inRange<T>(refK.z0, z0Min12, z0Max12) &&
            inRange<T>(refK.phiSlope, psMin12, psMax12))) {
        return;
      }
      refKUndecided = false;
    }
  }

  const int SP2 = M2_SP[k];
  const int m23SrcModule = doubletModule1[module23];
  const int srcWorkI =
      doubletSrcOffsets[module23] + (SP2 - hitIndices[m23SrcModule]);
  const int begin = edgeSumPerSrcHit[srcWorkI];
  const int end = edgeSumPerSrcHit[srcWorkI + 1];
  if (begin == end) {
    return;
  }

  // Load the per-triplet cuts once instead of on every loop iteration
  const T z0Min = MD23_z0_min[tripletIndex];
  const T z0Max = MD23_z0_max[tripletIndex];
  const T psMin = MD23_phi_slope_min[tripletIndex];
  const T psMax = MD23_phi_slope_max[tripletIndex];
  const T detaMin = MD23_deta_min[tripletIndex];
  const T detaMax = MD23_deta_max[tripletIndex];
  const T dphiMin = MD23_dphi_min[tripletIndex];
  const T dphiMax = MD23_dphi_max[tripletIndex];
  const T dydxMin = diffDydxMin[tripletIndex];
  const T dydxMax = diffDydxMax[tripletIndex];
  const T dzdrMin = diffDzdrMin[tripletIndex];
  const T dzdrMax = diffDzdrMax[tripletIndex];

  bool anyAccepted = false;

  for (int l = begin; l < end; ++l) {
    // ---- cuts on edge l (module doublet M2-M3) ----
    const float4 gl = geo[l];
    if (!(inRange<T>(gl.z, detaMin, detaMax) &&
          inRange<T>(gl.w, dphiMin, dphiMax))) {
      continue;
    }

    const float4 sl = edgeSlope[l];
    const int z0Ok = inRangeFiltered<T>(gl.x, sl.z, z0Min, z0Max);
    if (z0Ok == 0) {
      continue;
    }
    const int psOk =
        inRangeFiltered<T>(gl.y, kFilterErrScale * fabs(gl.y), psMin, psMax);
    if (psOk == 0) {
      continue;
    }

    // ---- cuts on the slope differences between k and l ----
    const int dydxOk = inRangeFiltered<T>(
        sk.x - sl.x, kFilterErrScale * (fabs(sk.x) + fabs(sl.x)), dydxMin,
        dydxMax);
    if (dydxOk == 0) {
      continue;
    }
    const int dzdrOk = inRangeFiltered<T>(
        sk.y - sl.y, kFilterErrScale * (fabs(sk.y) + fabs(sl.y)), dzdrMin,
        dzdrMax);
    if (dzdrOk == 0) {
      continue;
    }

    // ---- undecided cuts are resolved with the reference geometry ----
    if (refKUndecided || z0Ok < 0 || psOk < 0 || dydxOk < 0 || dzdrOk < 0) {
      if constexpr (kDeferFallback) {
        const int slot = atomicAdd(fallbackCount, 1);
        if (slot < fallbackCapacity) {
          fallbackPairs[slot] = TripletFallbackPair{tripletIndex, k, l};
        }
        continue;
      } else {
        const EdgeGeoReference<T> refL = edgeGeoReference<T>(
            M1_SP[l], M2_SP[l], R, z, x, y, phi, pi, epsilon);
        if (z0Ok < 0 && !inRange<T>(refL.z0, z0Min, z0Max)) {
          continue;
        }
        if (psOk < 0 && !inRange<T>(refL.phiSlope, psMin, psMax)) {
          continue;
        }
        if (dydxOk < 0 || dzdrOk < 0) {
          if (!haveRefK) {
            refK = edgeGeoReference<T>(M1_SP[k], M2_SP[k], R, z, x, y, phi, pi,
                                       epsilon);
            haveRefK = true;
          }
          if (dydxOk < 0 && !inRange<T>(static_cast<T>(refK.dydx - refL.dydx),
                                        dydxMin, dydxMax)) {
            continue;
          }
          if (dzdrOk < 0 && !inRange<T>(static_cast<T>(refK.dzdr - refL.dzdr),
                                        dzdrMin, dzdrMax)) {
            continue;
          }
        }
      }
    }

    edgeTag[l] = true;
    anyAccepted = true;
  }

  if (anyAccepted) {
    edgeTag[k] = true;
  }
}

// Resolves the pairs queued by tripletPairCutsFusedM23<T, true> with the
// reference geometry. Grid-stride loop, the queue length is only known on the
// device.
template <typename T>
__global__ void tripletPairCutsFallback(
    bool *edgeTag, const TripletFallbackPair *__restrict__ fallbackPairs,
    const int *__restrict__ fallbackCount, int fallbackCapacity,
    const float4 *__restrict__ geo, const T *__restrict__ MD12_z0_min,
    const T *__restrict__ MD12_phi_slope_min,
    const T *__restrict__ MD12_deta_min, const T *__restrict__ MD12_dphi_min,
    const T *__restrict__ MD12_z0_max, const T *__restrict__ MD12_phi_slope_max,
    const T *__restrict__ MD12_deta_max, const T *__restrict__ MD12_dphi_max,
    const T *__restrict__ MD23_z0_min, const T *__restrict__ MD23_phi_slope_min,
    const T *__restrict__ MD23_deta_min, const T *__restrict__ MD23_dphi_min,
    const T *__restrict__ MD23_z0_max, const T *__restrict__ MD23_phi_slope_max,
    const T *__restrict__ MD23_deta_max, const T *__restrict__ MD23_dphi_max,
    const T *__restrict__ diffDydxMin, const T *__restrict__ diffDydxMax,
    const T *__restrict__ diffDzdrMin, const T *__restrict__ diffDzdrMax,
    const int *__restrict__ M1_SP, const int *__restrict__ M2_SP,
    const T *__restrict__ R, const T *__restrict__ z, const T *__restrict__ x,
    const T *__restrict__ y, const T *__restrict__ phi, T pi, T epsilon) {
  const int n = min(*fallbackCount, fallbackCapacity);
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += gridDim.x * blockDim.x) {
    const TripletFallbackPair p = fallbackPairs[i];
    const int t = p.tripletIndex;
    const float4 gk = geo[p.k];
    const float4 gl = geo[p.l];
    const EdgeGeoReference<T> refK = edgeGeoReference<T>(
        M1_SP[p.k], M2_SP[p.k], R, z, x, y, phi, pi, epsilon);
    const EdgeGeoReference<T> refL = edgeGeoReference<T>(
        M1_SP[p.l], M2_SP[p.l], R, z, x, y, phi, pi, epsilon);

    const bool pass = inRange<T>(gk.z, MD12_deta_min[t], MD12_deta_max[t]) &&
                      inRange<T>(gk.w, MD12_dphi_min[t], MD12_dphi_max[t]) &&
                      inRange<T>(refK.z0, MD12_z0_min[t], MD12_z0_max[t]) &&
                      inRange<T>(refK.phiSlope, MD12_phi_slope_min[t],
                                 MD12_phi_slope_max[t]) &&
                      inRange<T>(gl.z, MD23_deta_min[t], MD23_deta_max[t]) &&
                      inRange<T>(gl.w, MD23_dphi_min[t], MD23_dphi_max[t]) &&
                      inRange<T>(refL.z0, MD23_z0_min[t], MD23_z0_max[t]) &&
                      inRange<T>(refL.phiSlope, MD23_phi_slope_min[t],
                                 MD23_phi_slope_max[t]) &&
                      inRange<T>(static_cast<T>(refK.dydx - refL.dydx),
                                 diffDydxMin[t], diffDydxMax[t]) &&
                      inRange<T>(static_cast<T>(refK.dzdr - refL.dzdr),
                                 diffDzdrMin[t], diffDzdrMax[t]);
    if (pass) {
      edgeTag[p.k] = true;
      edgeTag[p.l] = true;
    }
  }
}

}  // namespace ActsPlugins::detail
