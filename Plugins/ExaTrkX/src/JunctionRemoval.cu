// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"
#include "Acts/Plugins/ExaTrkX/detail/JunctionRemoval.hpp"

#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <thrust/transform_scan.h>

namespace Acts::detail {

__global__ void findNumInOutEdge(std::size_t nEdges, const float *scores,
                                 const std::int64_t *srcNodes,
                                 const std::int64_t *dstNodes, int *numInEdges,
                                 int *numOutEdges) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nEdges) {
    return;
  }

  int srcNode = srcNodes[i];
  int dstNode = dstNodes[i];

  atomicAdd(&numInEdges[dstNode], 1);
  atomicAdd(&numOutEdges[srcNode], 1);
}

__global__ void keepOnlyJunctions(std::size_t nNodes, int *numInEdges,
                                  int *numOutEdges) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nNodes) {
    return;
  }

  // Set the number of in and out edges to 0 for non-junction cases
  if (numInEdges[i] <= 1) {
    numInEdges[i] = 0;
  }
  if (numOutEdges[i] <= 1) {
    numOutEdges[i] = 0;
  }
}

__global__ void fillJunctionEdges(std::size_t nEdges,
                                  const std::int64_t *edgeNodes,
                                  const int *numEdgesPrefixSum,
                                  int *junctionEdges, int *junctionEdgeOffset) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nEdges) {
    return;
  }

  int node = edgeNodes[i];
  int base = numEdgesPrefixSum[node];
  int numEdgesNode = numEdgesPrefixSum[node + 1] - base;

  if (numEdgesNode == 1) {
    printf("[t%d, b%d] WARNING: node %d is not a junction\n", threadIdx.x,
           blockIdx.x, node);
    return;
  }

  if (numEdgesNode != 0) {
    int offset = atomicAdd(&junctionEdgeOffset[node], 1);
    if (offset >= numEdgesNode) {
      printf("[t%d, b%d] WARNING: inconsistent edge count for node %d\n",
             threadIdx.x, blockIdx.x, node);
      return;
    }
    junctionEdges[base + offset] = i;
  }
}

__global__ void fillEdgeMask(std::size_t nNodes, const float *scores,
                             const int *numEdgesPrefixSum,
                             const int *junctionEdges,
                             bool *edgesToRemoveMask) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nNodes) {
    return;
  }

  // Get the bse and number of edges for the current node
  int base = numEdgesPrefixSum[i];
  int numEdgesNode = numEdgesPrefixSum[i + 1] - base;

  // Find the edge with the maximum score
  float maxScore = 0.0f;
  int edgeIdMaxScore = -1;
  for (int j = base; j < base + numEdgesNode; ++j) {
    int edgeId = junctionEdges[j];
    float score = scores[edgeId];
    if (score > maxScore) {
      maxScore = score;
      edgeIdMaxScore = edgeId;
    }
  }

  // Mark all edges except the one with the maximum score for removal
  for (int j = base; j < base + numEdgesNode; ++j) {
    int edgeId = junctionEdges[j];
    if (edgeId != edgeIdMaxScore) {
      edgesToRemoveMask[edgeId] = true;
    }
  }
}

struct CastBoolToIntNot {
  int __device__ operator()(bool b) { return static_cast<int>(!b); }
};

struct LogicalNotPredicate {
  bool __device__ operator()(bool b) { return !b; }
};

std::pair<std::int64_t *, std::size_t> junctionRemovalCuda(
    std::size_t nEdges, std::size_t nNodes, const float *scores,
    const std::int64_t *srcNodes, const std::int64_t *dstNodes,
    cudaStream_t stream) {
  // Allocate device memory for the number of in and out edges
  int *numInEdges{}, *numOutEdges{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&numInEdges, (nNodes + 1) * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&numOutEdges, (nNodes + 1) * sizeof(int), stream));

  // Initialize the number of in and out edges to 0
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(numInEdges, 0, (nNodes + 1) * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(numOutEdges, 0, (nNodes + 1) * sizeof(int), stream));

  // Launch the kernel to find the number of in and out edges
  const dim3 blockSize = 512;
  const dim3 gridSizeEdges = (nEdges + blockSize.x - 1) / blockSize.x;
  findNumInOutEdge<<<gridSizeEdges, blockSize, 0, stream>>>(
      nEdges, scores, srcNodes, dstNodes, numInEdges, numOutEdges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Launch the kernel to keep only junctions
  const dim3 gridSizeNodes = (nNodes + blockSize.x - 1) / blockSize.x;
  keepOnlyJunctions<<<gridSizeNodes, blockSize, 0, stream>>>(nNodes, numInEdges,
                                                             numOutEdges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Perform prefix sum on the number of in and out edges
  thrust::exclusive_scan(thrust::device.on(stream), numInEdges,
                         numInEdges + nNodes + 1, numInEdges);
  thrust::exclusive_scan(thrust::device.on(stream), numOutEdges,
                         numOutEdges + nNodes + 1, numOutEdges);

  // Find the total number of in and out edges involved in junctions
  int numJunctionInEdges{}, numJunctionOutEdges{};
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&numJunctionInEdges, &numInEdges[nNodes],
                                  sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&numJunctionOutEdges, &numOutEdges[nNodes],
                                  sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  // Allocate device memory to store the edge indices for in and out edges
  int *junctionInEdges{}, *junctionOutEdges{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&junctionInEdges,
                                  numJunctionInEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&junctionOutEdges,
                                  numJunctionOutEdges * sizeof(int), stream));

  // Allocate device memory for the running index of the in and out edges per
  // node
  int *junctionInEdgeOffset{}, *junctionOutEdgeOffset{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&junctionInEdgeOffset, nNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&junctionOutEdgeOffset, nNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(junctionInEdgeOffset, 0, nNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(junctionOutEdgeOffset, 0, nNodes * sizeof(int), stream));

  // Fill the junction edges for in and out edges
  fillJunctionEdges<<<gridSizeEdges, blockSize, 0, stream>>>(
      nEdges, srcNodes, numOutEdges, junctionOutEdges, junctionOutEdgeOffset);
  ACTS_CUDA_CHECK(cudaGetLastError());
  fillJunctionEdges<<<gridSizeEdges, blockSize, 0, stream>>>(
      nEdges, dstNodes, numInEdges, junctionInEdges, junctionInEdgeOffset);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Allocate device memory for the edge mask
  bool *edgesToRemoveMask{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&edgesToRemoveMask, nEdges * sizeof(bool), stream));
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(edgesToRemoveMask, 0, nEdges * sizeof(bool), stream));

  // Fill the edge mask with the edges to be removed
  fillEdgeMask<<<gridSizeNodes, blockSize, 0, stream>>>(
      nNodes, scores, numInEdges, junctionInEdges, edgesToRemoveMask);
  ACTS_CUDA_CHECK(cudaGetLastError());
  fillEdgeMask<<<gridSizeNodes, blockSize, 0, stream>>>(
      nNodes, scores, numOutEdges, junctionOutEdges, edgesToRemoveMask);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Free the device memory
  ACTS_CUDA_CHECK(cudaFreeAsync(numInEdges, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(numOutEdges, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(junctionInEdges, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(junctionOutEdges, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(junctionInEdgeOffset, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(junctionOutEdgeOffset, stream));

  // Compactify the edges based on the edge mask
  int nEdgesToRemove =
      thrust::count(thrust::device.on(stream), edgesToRemoveMask,
                    edgesToRemoveMask + nEdges, true);
  int nEdgesAfter = nEdges - nEdgesToRemove;
  // Allocate memory for the new srcNodes and dstNodes arrays
  std::int64_t *newSrcNodes{};
  ACTS_CUDA_CHECK(cudaMallocAsync(
      &newSrcNodes, 2 * nEdgesAfter * sizeof(std::int64_t), stream));
  std::int64_t *newDstNodes = newSrcNodes + nEdgesAfter;

  // Compactify the srcNodes and dstNodes arrays based on the edge mask
  thrust::copy_if(thrust::device.on(stream), srcNodes, srcNodes + nEdges,
                  edgesToRemoveMask, newSrcNodes, LogicalNotPredicate{});
  thrust::copy_if(thrust::device.on(stream), dstNodes, dstNodes + nEdges,
                  edgesToRemoveMask, newDstNodes, LogicalNotPredicate{});

  // Free the device memory for the edge mask
  ACTS_CUDA_CHECK(cudaFreeAsync(edgesToRemoveMask, stream));

  // Synchronize the stream
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  return std::make_pair(newSrcNodes, static_cast<std::size_t>(nEdgesAfter));
}

}  // namespace Acts::detail
