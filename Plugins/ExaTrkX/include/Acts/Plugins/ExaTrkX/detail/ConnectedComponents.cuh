// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"

#include <cstdint>

#include <thrust/execution_policy.h>
#include <thrust/scan.h>

namespace Acts::detail {

/// Implementation of the FastSV algorithm as shown in
/// https://arxiv.org/abs/1910.05971

/// Hooking step of the FastSV algorithm
template <typename TEdge, typename TLabel>
__device__ void hookEdgesImpl(std::size_t i, const TEdge *sourceEdges,
                              const TEdge *targetEdges, const TLabel *labels,
                              TLabel *labelsNext, bool &changed) {
  auto u = sourceEdges[i];
  auto v = targetEdges[i];

  if (labels[u] == labels[labels[u]] && labels[v] < labels[u]) {
    labelsNext[labels[u]] = labels[v];
    changed = true;
    //printf("Edge (%i,%i): set labelsNext[%i] = labels[%i] = %i\n", u, v, labels[u], v, labels[v]);
  } else if (labels[v] == labels[labels[v]] && labels[u] < labels[v]) {
    labelsNext[labels[v]] = labels[u];
    changed = true;
    //printf("Edge (%i,%i): set labelsNext[%i] = labels[%i] = %i\n", u, v, labels[v], u, labels[u]);
  } else {
    //printf("Edge (%i,%i): no action\n", u, v);
  }
}

/// Shortcutting step of the FastSV algorithm
template <typename TEdge, typename TLabel>
__device__ void shortcutImpl(std::size_t i, const TEdge *sourceEdges,
                             const TEdge *targetEdges, const TLabel *labels,
                             TLabel *labelsNext, bool &changed) {
  if (labels[i] != labels[labels[i]]) {
    labelsNext[i] = labels[labels[i]];
    //printf("Vertex %i: labelsNext[%i] = labels[%i] = %i\n", i, i, labels[i], labels[labels[i]]);
    changed = true;
  }
}

/// Implementation of the FastSV algorithm in a single kernel
/// NOTE: This can only run in one block due to synchronization
template <typename TEdge, typename TLabel>
__global__ void labelConnectedComponents(std::size_t numEdges,
                                         const TEdge *sourceEdges,
                                         const TEdge *targetEdges,
                                         std::size_t numNodes, TLabel *labels,
                                         TLabel *labelsNext) {
  // Currently this kernel works only with 1 block
  assert(gridDim.x == 1 && gridDim.y == 1 && gridDim.z == 1);

  for (std::size_t i = threadIdx.x; i < numNodes; i += blockDim.x) {
    labels[i] = i;
    labelsNext[i] = i;
  }

  __syncthreads();
  bool changed = false;

  do {
    changed = false;

    //printf("Iteration %i\n", n);

    // Tree hooking for each edge;
    for (std::size_t i = threadIdx.x; i < numEdges; i += blockDim.x) {
      hookEdgesImpl(i, sourceEdges, targetEdges, labels, labelsNext, changed);
    }
    __syncthreads();

    for (std::size_t i = threadIdx.x; i < numNodes; i += blockDim.x) {
      labels[i] = labelsNext[i];
    }

    /*if(threadIdx.x == 0 ) {
      for(int i=0; i<numNodes; ++i) {
        printf("Vertex %i - label %i\n", i, labels[i]);
      }
    }*/

    // Shortcutting
    for (std::size_t i = threadIdx.x; i < numNodes; i += blockDim.x) {
      shortcutImpl(i, sourceEdges, targetEdges, labels, labelsNext, changed);
    }

    for (std::size_t i = threadIdx.x; i < numNodes; i += blockDim.x) {
      labels[i] = labelsNext[i];
    }

    /*if(threadIdx.x == 0 ) {
      for(int i=0; i<numNodes; ++i) {
        printf("Vertex after Shortcutting %i - label %i\n", i, labels[i]);
      }
    }*/

  } while (__syncthreads_or(changed));
}

/// Hooking-kernel for implementing the FastSV loop on the host
template <typename TEdge, typename TLabel>
__global__ void hookEdges(std::size_t numEdges, const TEdge *sourceEdges,
                          const TEdge *targetEdges, const TLabel *labels,
                          TLabel *labelsNext, int *globalChanged) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  bool changed = false;
  if (i < numEdges) {
    hookEdgesImpl(i, sourceEdges, targetEdges, labels, labelsNext, changed);
  }

  if (__syncthreads_or(changed) && threadIdx.x == 0) {
    *globalChanged = true;
  }
}

/// Shortcutting-kernel for implementing the FastSV loop on the host
template <typename TEdge, typename TLabel>
__global__ void shortcut(std::size_t numNodes, const TEdge *sourceEdges,
                         const TEdge *targetEdges, const TLabel *labels,
                         TLabel *labelsNext, int *globalChanged) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  bool changed = false;
  if (i < numNodes) {
    shortcutImpl(i, sourceEdges, targetEdges, labels, labelsNext, changed);
  }

  if (__syncthreads_or(changed) && threadIdx.x == 0) {
    *globalChanged = true;
  }
}

template <typename T>
__global__ void makeLabelMask(std::size_t nLabels, const T *labels,
                              T *labelMask) {
  std::size_t i = threadIdx.x + blockDim.x * blockIdx.x;

  if (i >= nLabels) {
    return;
  }

  labelMask[labels[i]] = 1;
}

template <typename T>
__global__ void mapEdgeLabels(std::size_t nLabels, T *labels,
                              const T *mapping) {
  std::size_t i = threadIdx.x + blockDim.x * blockIdx.x;

  if (i >= nLabels) {
    return;
  }

  labels[i] = mapping[labels[i]];
}

template <typename TEdges, typename TLabel>
TLabel connectedComponentsCuda(std::size_t nEdges, const TEdges *sourceEdges,
                               const TEdges *targetEdges, std::size_t nNodes,
                               TLabel *labels, cudaStream_t stream,
                               bool useOneCudaBlock = true) {
  TLabel *tmpMemory = nullptr;
  ACTS_CUDA_CHECK(cudaMallocAsync(&tmpMemory, nNodes * sizeof(TLabel), stream));

  const dim3 blockDim = 1024;

  if (useOneCudaBlock) {
    // Make synchronization in one block, to avoid that inter-block sync is
    // necessary
    labelConnectedComponents<<<1, blockDim, 1, stream>>>(
        nEdges, sourceEdges, targetEdges, nNodes, labels, tmpMemory);
    ACTS_CUDA_CHECK(cudaGetLastError());
  } else {
    int changed = false;
    int *cudaChanged = nullptr;
    ACTS_CUDA_CHECK(cudaMallocAsync(&cudaChanged, sizeof(int), stream));

    const dim3 gridDimNodes = (nNodes + blockDim.x - 1) / blockDim.x;
    const dim3 gridDimEdges = (nEdges + blockDim.x - 1) / blockDim.x;

    detail::iota<<<gridDimNodes, blockDim, 0, stream>>>(nNodes, labels);
    ACTS_CUDA_CHECK(cudaGetLastError());
    detail::iota<<<gridDimNodes, blockDim, 0, stream>>>(nNodes, tmpMemory);
    ACTS_CUDA_CHECK(cudaGetLastError());

    do {
      ACTS_CUDA_CHECK(cudaMemsetAsync(cudaChanged, 0, sizeof(int), stream));

      // Hooking
      hookEdges<<<gridDimEdges, blockDim, 0, stream>>>(
          nEdges, sourceEdges, targetEdges, labels, tmpMemory, cudaChanged);
      ACTS_CUDA_CHECK(cudaGetLastError());
      ACTS_CUDA_CHECK(cudaMemcpyAsync(labels, tmpMemory,
                                      nNodes * sizeof(TLabel),
                                      cudaMemcpyDeviceToDevice, stream));

      // Shortcutting
      shortcut<<<gridDimNodes, blockDim, 0, stream>>>(
          nNodes, sourceEdges, targetEdges, labels, tmpMemory, cudaChanged);
      ACTS_CUDA_CHECK(cudaGetLastError());
      ACTS_CUDA_CHECK(cudaMemcpyAsync(labels, tmpMemory,
                                      nNodes * sizeof(TLabel),
                                      cudaMemcpyDeviceToDevice, stream));

      ACTS_CUDA_CHECK(cudaMemcpyAsync(&changed, cudaChanged, sizeof(int),
                                      cudaMemcpyDeviceToHost, stream));
      ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    } while (changed);

    ACTS_CUDA_CHECK(cudaFreeAsync(cudaChanged, stream));
  }

  // Assume we have the following components:
  // 0 3 5 3 0 0

  // Fill a mask which labels survived the connected components algorithm
  // 0 1 2 3 4 5
  // 1 0 0 1 0 1
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(tmpMemory, 0, nNodes * sizeof(TLabel), stream));
  dim3 gridDim = (nNodes + blockDim.x - 1) / blockDim.x;
  makeLabelMask<<<gridDim, blockDim, 0, stream>>>(nNodes, labels, tmpMemory);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Exclusive prefix sum on the label mask
  // 0 1 2 3 4 5
  // 0 1 1 1 2 2
  thrust::exclusive_scan(thrust::device.on(stream), tmpMemory,
                         tmpMemory + nNodes, tmpMemory);

  // Remap edge labels with values in prefix sum
  // 0 -> 0, 3 -> 1, 5 -> 2
  mapEdgeLabels<<<gridDim, blockDim, 0, stream>>>(nNodes, labels, tmpMemory);
  ACTS_CUDA_CHECK(cudaGetLastError());

  TLabel nLabels;
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&nLabels, &tmpMemory[nNodes - 1],
                                  sizeof(TLabel), cudaMemcpyDeviceToHost,
                                  stream));

  ACTS_CUDA_CHECK(cudaFreeAsync(tmpMemory, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  return nLabels;
}

}  // namespace Acts::detail
