// This file is part of the ACTS project
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/DWalkTrackBuilding.hpp"

#include "ActsPlugins/Gnn/detail/ConnectedComponents.cuh"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <cuda_runtime_api.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kBlockSize = 256;
constexpr std::size_t kSmallResidualEdgeThreshold = 512;

struct IsNonZero {
  __host__ __device__ bool operator()(unsigned char value) const {
    return value != 0;
  }
};

struct Edge {
  int src = 0;
  int dst = 0;
};

struct DeviceOrientedEdges {
  int *src = nullptr;
  int *dst = nullptr;
  float *score = nullptr;
  unsigned char *valid = nullptr;
  unsigned char *activeNodes = nullptr;
  std::size_t numEdges = 0;
  std::size_t numNodes = 0;
};

struct DeviceCompactEdges {
  int *src = nullptr;
  int *dst = nullptr;
  float *score = nullptr;
  std::size_t numEdges = 0;
  std::size_t numNodes = 0;
};

struct DeviceCsrGraph {
  int *rowPtr = nullptr;
  int *colIdx = nullptr;
  float *edgeWeight = nullptr;
  int *incomingRowPtr = nullptr;
  int *incomingColIdx = nullptr;
  std::size_t numNodes = 0;
  std::size_t numEdges = 0;
};

struct DpCudaState {
  float *bestScore = nullptr;
  int *bestChild = nullptr;
  unsigned char *sourceMask = nullptr;
  std::size_t numNodes = 0;
};

void freeDeviceOrientedEdges(DeviceOrientedEdges &graph, cudaStream_t stream) {
  if (graph.src != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.src, stream));
  }
  if (graph.dst != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.dst, stream));
  }
  if (graph.score != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.score, stream));
  }
  if (graph.valid != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.valid, stream));
  }
  if (graph.activeNodes != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.activeNodes, stream));
  }
  graph = {};
}

void freeDeviceCompactEdges(DeviceCompactEdges &graph, cudaStream_t stream) {
  if (graph.src != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.src, stream));
  }
  if (graph.dst != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.dst, stream));
  }
  if (graph.score != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.score, stream));
  }
  graph = {};
}

void freeDeviceCsrGraph(DeviceCsrGraph &graph, cudaStream_t stream) {
  if (graph.rowPtr != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.rowPtr, stream));
  }
  if (graph.colIdx != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.colIdx, stream));
  }
  if (graph.edgeWeight != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.edgeWeight, stream));
  }
  if (graph.incomingRowPtr != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.incomingRowPtr, stream));
  }
  if (graph.incomingColIdx != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(graph.incomingColIdx, stream));
  }
  graph = {};
}

void freeDpCudaState(DpCudaState &state, cudaStream_t stream) {
  if (state.bestScore != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(state.bestScore, stream));
  }
  if (state.bestChild != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(state.bestChild, stream));
  }
  if (state.sourceMask != nullptr) {
    ACTS_CUDA_CHECK(cudaFreeAsync(state.sourceMask, stream));
  }
  state = {};
}

std::vector<int> orderedSimpleComponentNodes(const std::vector<int> &nodes,
                                             int component,
                                             const std::vector<int> &labels,
                                             const std::vector<int> &inDegree,
                                             const std::vector<int> &nextNode) {
  if (nodes.empty()) {
    return {};
  }

  int start = nodes.front();
  for (int node : nodes) {
    if (inDegree.at(node) == 0) {
      start = node;
      break;
    }
  }

  std::vector<int> ordered;
  ordered.reserve(nodes.size());
  int node = start;
  while (node >= 0 && labels.at(node) == component &&
         ordered.size() < nodes.size()) {
    ordered.push_back(node);
    node = nextNode.at(node);
  }

  if (ordered.size() != nodes.size()) {
    ordered = nodes;
    std::sort(ordered.begin(), ordered.end());
  }
  return ordered;
}

float minRootScore(const std::string &pathMetric) {
  if (pathMetric == "score_weighted_length") {
    return 0.0F;
  }
  if (pathMetric == "length") {
    return 2.0F;
  }
  throw std::invalid_argument(
      "DWalkTrackBuilding pathMetric must be 'score_weighted_length' or "
      "'length'");
}

class DisjointSet {
 public:
  explicit DisjointSet(std::size_t size) : m_parent(size), m_rank(size, 0) {
    for (std::size_t i = 0; i < size; ++i) {
      m_parent[i] = static_cast<int>(i);
    }
  }

  int find(int value) {
    int parent = m_parent.at(value);
    if (parent != value) {
      parent = find(parent);
      m_parent.at(value) = parent;
    }
    return parent;
  }

  void unite(int lhs, int rhs) {
    int lhsRoot = find(lhs);
    int rhsRoot = find(rhs);
    if (lhsRoot == rhsRoot) {
      return;
    }
    if (m_rank.at(lhsRoot) < m_rank.at(rhsRoot)) {
      std::swap(lhsRoot, rhsRoot);
    }
    m_parent.at(rhsRoot) = lhsRoot;
    if (m_rank.at(lhsRoot) == m_rank.at(rhsRoot)) {
      ++m_rank.at(lhsRoot);
    }
  }

 private:
  std::vector<int> m_parent;
  std::vector<int> m_rank;
};

__device__ void atomicMaxFloat(float *address, float value) {
  int *addressAsInt = reinterpret_cast<int *>(address);
  int old = *addressAsInt;
  int assumed;
  do {
    assumed = old;
    old = atomicCAS(addressAsInt, assumed,
                    __float_as_int(fmaxf(value, __int_as_float(assumed))));
  } while (assumed != old);
}

__global__ void orientEdgesKernel(std::size_t numEdges, std::size_t numNodes,
                                  std::size_t numFeatures,
                                  std::size_t radialFeatureIndex,
                                  const std::int64_t *edgeIndex,
                                  const float *scores,
                                  const float *nodeFeatures, int *srcOut,
                                  int *dstOut, float *scoreOut,
                                  unsigned char *validEdge) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }

  int src = static_cast<int>(edgeIndex[edge]);
  int dst = static_cast<int>(edgeIndex[numEdges + edge]);
  if (src < 0 || dst < 0 || src >= static_cast<int>(numNodes) ||
      dst >= static_cast<int>(numNodes) || src == dst) {
    srcOut[edge] = 0;
    dstOut[edge] = 0;
    scoreOut[edge] = 0.0F;
    validEdge[edge] = 0;
    return;
  }

  float srcR = nodeFeatures[src * numFeatures + radialFeatureIndex];
  float dstR = nodeFeatures[dst * numFeatures + radialFeatureIndex];
  if (srcR > dstR || (srcR == dstR && src > dst)) {
    int tmp = src;
    src = dst;
    dst = tmp;
  }

  srcOut[edge] = src;
  dstOut[edge] = dst;
  scoreOut[edge] = scores[edge];
  validEdge[edge] = 1;
}

__global__ void markActiveNodesKernel(std::size_t numEdges, const int *src,
                                      const int *dst,
                                      const unsigned char *validEdge,
                                      unsigned char *activeNodes) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges || validEdge[edge] == 0) {
    return;
  }
  activeNodes[src[edge]] = 1;
  activeNodes[dst[edge]] = 1;
}

__global__ void computeComponentDegreeKernel(
    std::size_t numEdges, const int *src, const int *dst,
    const unsigned char *validEdge, int *inDegree, int *outDegree) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges || validEdge[edge] == 0) {
    return;
  }

  int s = src[edge];
  int d = dst[edge];
  atomicAdd(outDegree + s, 1);
  atomicAdd(inDegree + d, 1);
}

__global__ void countActiveComponentNodesKernel(
    std::size_t numNodes, const unsigned char *activeNodes, const int *labels,
    int *componentSizes) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || activeNodes[node] == 0) {
    return;
  }
  atomicAdd(componentSizes + labels[node], 1);
}

__global__ void markBadComponentsKernel(
    std::size_t numNodes, const unsigned char *activeNodes, const int *labels,
    const int *inDegree, const int *outDegree, int *badComponents) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || activeNodes[node] == 0) {
    return;
  }

  if (max(inDegree[node], outDegree[node]) > 1) {
    badComponents[labels[node]] = 1;
  }
}

__global__ void buildInitialComponentMasksKernel(
    std::size_t numNodes, const unsigned char *activeNodes, const int *labels,
    const int *componentSizes, const int *badComponents,
    int minCandidateSize, unsigned char *simpleNodeMask,
    unsigned char *complexNodeMask) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || activeNodes[node] == 0) {
    if (node < numNodes) {
      simpleNodeMask[node] = 0;
      complexNodeMask[node] = 0;
    }
    return;
  }

  int label = labels[node];
  bool isLarge = componentSizes[label] >= minCandidateSize;
  bool isSimple = isLarge && badComponents[label] == 0;
  simpleNodeMask[node] = isSimple ? 1 : 0;
  complexNodeMask[node] = (isLarge && !isSimple) ? 1 : 0;
}

__global__ void maskEdgesByActiveNodesKernel(std::size_t numEdges,
                                             const int *src, const int *dst,
                                             const unsigned char *activeNodes,
                                             unsigned char *keepMask) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  keepMask[edge] = activeNodes[src[edge]] && activeNodes[dst[edge]] ? 1 : 0;
}

__global__ void deactivateSelectedNodesKernel(std::size_t numSelected,
                                              const int *selectedNodes,
                                              unsigned char *activeNodes) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= numSelected) {
    return;
  }
  int node = selectedNodes[i];
  if (node >= 0) {
    activeNodes[node] = 0;
  }
}

__global__ void fillCountsKernel(std::size_t numEdges, const int *nodes,
                                 int *counts) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  atomicAdd(counts + nodes[edge], 1);
}

__global__ void scatterSortedOutgoingKernel(
    std::size_t numEdges, const int *sortedSrc, const int *sortedDst,
    const float *sortedScore, const int *rowPtr, int *cursor, int *colIdx,
    float *edgeWeight, bool useLengthMetric) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  int src = sortedSrc[edge];
  int slot = atomicAdd(cursor + src, 1);
  colIdx[rowPtr[src] + slot] = sortedDst[edge];
  edgeWeight[rowPtr[src] + slot] = useLengthMetric ? 1.0F : sortedScore[edge];
}

__global__ void scatterSortedIncomingKernel(std::size_t numEdges,
                                            const int *sortedDst,
                                            const int *sortedSrc,
                                            const int *incomingRowPtr,
                                            int *cursor,
                                            int *incomingColIdx) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  int dst = sortedDst[edge];
  int slot = atomicAdd(cursor + dst, 1);
  incomingColIdx[incomingRowPtr[dst] + slot] = sortedSrc[edge];
}

__global__ void initFloatKernel(std::size_t size, float value, float *array) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    array[i] = value;
  }
}

__global__ void maxAddBestScoresKernel(std::size_t numEdges, const int *src,
                                       const int *dst, const float *score,
                                       const unsigned char *validEdge,
                                       const unsigned char *activeNodes,
                                       float *bestOutScore,
                                       float *bestInScore) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  if (validEdge[edge] == 0) {
    return;
  }
  int s = src[edge];
  int d = dst[edge];
  if (!activeNodes[s] || !activeNodes[d]) {
    return;
  }
  atomicMaxFloat(bestOutScore + s, score[edge]);
  atomicMaxFloat(bestInScore + d, score[edge]);
}

__global__ void maxAddMaskKernel(std::size_t numEdges, const int *src,
                                 const int *dst, const float *score,
                                 const unsigned char *validEdge,
                                 const unsigned char *activeNodes,
                                 const float *bestOutScore,
                                 const float *bestInScore, float thMin,
                                 float thAdd, unsigned char *keepMask) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  if (validEdge[edge] == 0) {
    keepMask[edge] = 0;
    return;
  }

  int s = src[edge];
  int d = dst[edge];
  if (!activeNodes[s] || !activeNodes[d]) {
    keepMask[edge] = 0;
    return;
  }

  float edgeScore = score[edge];
  bool maskAdd = edgeScore > thAdd;
  bool outgoingKeep =
      bestOutScore[s] >= thMin && edgeScore == bestOutScore[s];
  bool incomingKeep = bestInScore[d] >= thMin && edgeScore == bestInScore[d];
  keepMask[edge] = ((outgoingKeep || maskAdd) && (incomingKeep || maskAdd)) ? 1
                                                                            : 0;
}

__global__ void initializeDpKernel(float *bestScore, int *bestChild,
                                  unsigned char *sourceMask, int *inDegree,
                                  int *remainingOutDegree,
                                  std::size_t numNodes) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes) {
    return;
  }
  bestScore[node] = 0.0F;
  bestChild[node] = -1;
  sourceMask[node] = 0;
  inDegree[node] = 0;
  remainingOutDegree[node] = 0;
}

__global__ void computeActiveDegreeKernel(const int *rowPtr, const int *colIdx,
                                          const unsigned char *activeNodes,
                                          int *inDegree,
                                          int *remainingOutDegree,
                                          std::size_t numNodes) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || !activeNodes[node]) {
    return;
  }

  int activeOutDegree = 0;
  for (int edge = rowPtr[node]; edge < rowPtr[node + 1]; ++edge) {
    int child = colIdx[edge];
    if (activeNodes[child]) {
      ++activeOutDegree;
      atomicAdd(inDegree + child, 1);
    }
  }
  remainingOutDegree[node] = activeOutDegree;
}

__global__ void initializeFrontierKernel(const unsigned char *activeNodes,
                                         const int *remainingOutDegree,
                                         int *frontier, int *frontierSize,
                                         std::size_t numNodes) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes) {
    return;
  }
  if (activeNodes[node] && remainingOutDegree[node] == 0) {
    int slot = atomicAdd(frontierSize, 1);
    frontier[slot] = static_cast<int>(node);
  }
}

__global__ void processFrontierKernel(const int *frontier, int frontierSize,
                                      const int *rowPtr, const int *colIdx,
                                      const float *edgeWeight,
                                      const unsigned char *activeNodes,
                                      const float *bestScore,
                                      float *frontierBestScore,
                                      int *frontierBestChild) {
  int frontierIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (frontierIndex >= frontierSize) {
    return;
  }

  int node = frontier[frontierIndex];
  float bestValue = 0.0F;
  int bestNext = -1;
  for (int edge = rowPtr[node]; edge < rowPtr[node + 1]; ++edge) {
    int child = colIdx[edge];
    if (!activeNodes[child]) {
      continue;
    }
    float candidate = edgeWeight[edge] + bestScore[child];
    if (candidate > bestValue) {
      bestValue = candidate;
      bestNext = child;
    }
  }

  frontierBestScore[frontierIndex] = bestValue;
  frontierBestChild[frontierIndex] = bestNext;
}

__global__ void finalizeFrontierKernel(const int *frontier, int frontierSize,
                                       const float *frontierBestScore,
                                       const int *frontierBestChild,
                                       float *bestScore, int *bestChild) {
  int frontierIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (frontierIndex >= frontierSize) {
    return;
  }

  int node = frontier[frontierIndex];
  bestScore[node] = frontierBestScore[frontierIndex];
  bestChild[node] = frontierBestChild[frontierIndex];
}

__global__ void enqueueParentFrontierKernel(
    const int *frontier, int frontierSize, const int *incomingRowPtr,
    const int *incomingColIdx, const unsigned char *activeNodes,
    int *remainingOutDegree, int *nextFrontier, int *nextFrontierSize) {
  int frontierIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (frontierIndex >= frontierSize) {
    return;
  }

  int node = frontier[frontierIndex];
  for (int edge = incomingRowPtr[node]; edge < incomingRowPtr[node + 1];
       ++edge) {
    int parent = incomingColIdx[edge];
    if (!activeNodes[parent]) {
      continue;
    }
    int oldValue = atomicSub(remainingOutDegree + parent, 1);
    if (oldValue == 1) {
      int slot = atomicAdd(nextFrontierSize, 1);
      nextFrontier[slot] = parent;
    }
  }
}

__global__ void buildSourceMaskKernel(const unsigned char *activeNodes,
                                      const int *inDegree,
                                      const float *bestScore,
                                      unsigned char *sourceMask,
                                      std::size_t numNodes) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes) {
    return;
  }
  sourceMask[node] =
      activeNodes[node] && inDegree[node] == 0 && bestScore[node] > 0.0F ? 1
                                                                         : 0;
}

__global__ void selectComponentScoresKernel(
    std::size_t numNodes, const int *componentLabels,
    const unsigned char *sourceMask, const float *bestScore,
    float *componentBestScore) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || sourceMask[node] == 0) {
    return;
  }
  int component = componentLabels[node];
  if (component < 0) {
    return;
  }
  atomicMaxFloat(componentBestScore + component, bestScore[node]);
}

__global__ void selectRootsKernel(std::size_t numNodes,
                                  const int *componentLabels,
                                  const unsigned char *sourceMask,
                                  const float *bestScore,
                                  const float *componentBestScore,
                                  float minRootScore, int *selectedRoots) {
  std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
  if (node >= numNodes || sourceMask[node] == 0) {
    return;
  }
  int component = componentLabels[node];
  if (component < 0 || componentBestScore[component] <= minRootScore) {
    return;
  }
  if (bestScore[node] == componentBestScore[component]) {
    atomicCAS(selectedRoots + component, -1, static_cast<int>(node));
  }
}

__global__ void traceSelectedPathsKernel(const int *bestChild,
                                         const int *selectedRoots,
                                         int *selectedTrackLabels,
                                         int *selectedNodes,
                                         int *selectedCount,
                                         std::size_t numNodes,
                                         std::size_t numPaths) {
  std::size_t pathIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (pathIndex >= numPaths) {
    return;
  }

  int node = selectedRoots[pathIndex];
  if (node < 0) {
    return;
  }
  std::size_t step = 0;
  while (node != -1 && step < numNodes) {
    int slot = atomicAdd(selectedCount, 1);
    selectedTrackLabels[slot] = static_cast<int>(pathIndex);
    selectedNodes[slot] = node;
    node = bestChild[node];
    ++step;
  }
}

__global__ void remapEdgesKernel(std::size_t numEdges, int *src, int *dst,
                                 const int *oldToNew) {
  std::size_t edge = blockIdx.x * blockDim.x + threadIdx.x;
  if (edge >= numEdges) {
    return;
  }
  src[edge] = oldToNew[src[edge]];
  dst[edge] = oldToNew[dst[edge]];
}

std::vector<Edge> createOrientedEdgesCuda(
    const ActsPlugins::Tensor<std::int64_t> &edgeTensor,
    const ActsPlugins::Tensor<float> &scoreTensor,
    const ActsPlugins::Tensor<float> &featureTensor,
    std::size_t radialFeatureIndex, cudaStream_t stream,
    DeviceOrientedEdges &deviceGraph, int **cudaInitialLabels = nullptr,
    int *initialNumComponents = nullptr,
    std::vector<int> *initialLabels = nullptr) {
  const auto numNodes = featureTensor.shape().at(0);
  const auto numFeatures = featureTensor.shape().at(1);
  const auto numEdges = edgeTensor.shape().at(1);

  deviceGraph.numEdges = numEdges;
  deviceGraph.numNodes = numNodes;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&deviceGraph.src, numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&deviceGraph.dst, numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&deviceGraph.score, numEdges * sizeof(float), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&deviceGraph.valid, numEdges * sizeof(unsigned char),
                      stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&deviceGraph.activeNodes,
                                  numNodes * sizeof(unsigned char), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(deviceGraph.activeNodes, 0,
                                  numNodes * sizeof(unsigned char), stream));

  const dim3 edgeGrid((numEdges + kBlockSize - 1) / kBlockSize);
  orientEdgesKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      numEdges, numNodes, numFeatures, radialFeatureIndex, edgeTensor.data(),
      scoreTensor.data(), featureTensor.data(), deviceGraph.src,
      deviceGraph.dst, deviceGraph.score, deviceGraph.valid);
  ACTS_CUDA_CHECK(cudaGetLastError());
  markActiveNodesKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      numEdges, deviceGraph.src, deviceGraph.dst, deviceGraph.valid,
      deviceGraph.activeNodes);
  ACTS_CUDA_CHECK(cudaGetLastError());

  if (cudaInitialLabels != nullptr || initialLabels != nullptr ||
      initialNumComponents != nullptr) {
    int *cudaLabels{};
    ACTS_CUDA_CHECK(
        cudaMallocAsync(&cudaLabels, numNodes * sizeof(int), stream));
    int numComponents = ActsPlugins::detail::connectedComponentsCuda(
        numEdges, deviceGraph.src, deviceGraph.dst, numNodes, cudaLabels,
        stream, false);
    if (initialNumComponents != nullptr) {
      *initialNumComponents = numComponents;
    }
    if (initialLabels != nullptr) {
      initialLabels->assign(numNodes, -1);
      ACTS_CUDA_CHECK(cudaMemcpyAsync(initialLabels->data(), cudaLabels,
                                      numNodes * sizeof(int),
                                      cudaMemcpyDeviceToHost, stream));
    }
    if (cudaInitialLabels != nullptr) {
      *cudaInitialLabels = cudaLabels;
    } else {
      ACTS_CUDA_CHECK(cudaFreeAsync(cudaLabels, stream));
    }
  }

  std::vector<int> src(numEdges);
  std::vector<int> dst(numEdges);
  std::vector<unsigned char> valid(numEdges);

  ACTS_CUDA_CHECK(cudaMemcpyAsync(src.data(), deviceGraph.src,
                                  numEdges * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(dst.data(), deviceGraph.dst,
                                  numEdges * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(valid.data(), deviceGraph.valid,
                                  numEdges * sizeof(unsigned char),
                                  cudaMemcpyDeviceToHost, stream));

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  std::vector<Edge> edges;
  edges.reserve(numEdges);
  for (std::size_t edge = 0; edge < numEdges; ++edge) {
    if (valid.at(edge) != 0) {
      edges.push_back({src.at(edge), dst.at(edge)});
    }
  }
  return edges;
}

DeviceCompactEdges compactDeviceEdges(const int *src, const int *dst,
                                      const float *score,
                                      const unsigned char *keepMask,
                                      std::size_t numEdges,
                                      std::size_t numNodes,
                                      cudaStream_t stream) {
  DeviceCompactEdges compact;
  compact.numNodes = numNodes;
  compact.numEdges = thrust::count(thrust::device.on(stream), keepMask,
                                   keepMask + numEdges, 1);
  if (compact.numEdges == 0) {
    return compact;
  }

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&compact.src, compact.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&compact.dst, compact.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&compact.score,
                                  compact.numEdges * sizeof(float), stream));

  thrust::copy_if(thrust::device.on(stream), src, src + numEdges, keepMask,
                  compact.src, IsNonZero{});
  thrust::copy_if(thrust::device.on(stream), dst, dst + numEdges, keepMask,
                  compact.dst, IsNonZero{});
  thrust::copy_if(thrust::device.on(stream), score, score + numEdges, keepMask,
                  compact.score, IsNonZero{});
  return compact;
}

DeviceCompactEdges maxAddCompactDeviceEdgesCuda(
    const DeviceOrientedEdges &deviceGraph, const unsigned char *activeNodes,
    float thMin, float thAdd, cudaStream_t stream) {
  if (deviceGraph.numEdges == 0) {
    return {};
  }

  float *cudaBestOut{};
  float *cudaBestIn{};
  unsigned char *cudaKeep{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaBestOut,
                                  deviceGraph.numNodes * sizeof(float),
                                  stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaBestIn,
                                  deviceGraph.numNodes * sizeof(float),
                                  stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaKeep,
                                  deviceGraph.numEdges * sizeof(unsigned char),
                                  stream));

  const dim3 nodeGrid((deviceGraph.numNodes + kBlockSize - 1) / kBlockSize);
  const dim3 edgeGrid((deviceGraph.numEdges + kBlockSize - 1) / kBlockSize);
  initFloatKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numNodes, -std::numeric_limits<float>::infinity(),
      cudaBestOut);
  initFloatKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numNodes, -std::numeric_limits<float>::infinity(),
      cudaBestIn);
  ACTS_CUDA_CHECK(cudaGetLastError());
  maxAddBestScoresKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numEdges, deviceGraph.src, deviceGraph.dst,
      deviceGraph.score, deviceGraph.valid, activeNodes, cudaBestOut,
      cudaBestIn);
  ACTS_CUDA_CHECK(cudaGetLastError());
  maxAddMaskKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numEdges, deviceGraph.src, deviceGraph.dst,
      deviceGraph.score, deviceGraph.valid, activeNodes, cudaBestOut,
      cudaBestIn, thMin, thAdd, cudaKeep);
  ACTS_CUDA_CHECK(cudaGetLastError());

  auto compact =
      compactDeviceEdges(deviceGraph.src, deviceGraph.dst, deviceGraph.score,
                         cudaKeep, deviceGraph.numEdges, deviceGraph.numNodes,
                         stream);

  ACTS_CUDA_CHECK(cudaFreeAsync(cudaBestOut, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaBestIn, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaKeep, stream));
  return compact;
}

DeviceCompactEdges compactActiveEdgesCuda(const DeviceCompactEdges &edges,
                                          const unsigned char *activeNodes,
                                          cudaStream_t stream) {
  if (edges.numEdges == 0) {
    return {};
  }

  unsigned char *cudaKeep{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaKeep,
                                  edges.numEdges * sizeof(unsigned char),
                                  stream));
  const dim3 edgeGrid((edges.numEdges + kBlockSize - 1) / kBlockSize);
  maskEdgesByActiveNodesKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      edges.numEdges, edges.src, edges.dst, activeNodes, cudaKeep);
  ACTS_CUDA_CHECK(cudaGetLastError());
  auto compact = compactDeviceEdges(edges.src, edges.dst, edges.score,
                                    cudaKeep, edges.numEdges, edges.numNodes,
                                    stream);
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaKeep, stream));
  return compact;
}

void classifyInitialComponentsCuda(
    const DeviceOrientedEdges &deviceGraph, const int *cudaLabels,
    int numComponents, int minCandidateSize, unsigned char **cudaSimpleNodeMask,
    unsigned char **cudaComplexNodeMask, cudaStream_t stream) {
  int *cudaInDegree{};
  int *cudaOutDegree{};
  int *cudaComponentSizes{};
  int *cudaBadComponents{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaInDegree,
                                  deviceGraph.numNodes * sizeof(int),
                                  stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaOutDegree,
                                  deviceGraph.numNodes * sizeof(int),
                                  stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaComponentSizes,
                                  numComponents * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaBadComponents,
                                  numComponents * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(cudaSimpleNodeMask,
                                  deviceGraph.numNodes * sizeof(unsigned char),
                                  stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(cudaComplexNodeMask,
                                  deviceGraph.numNodes * sizeof(unsigned char),
                                  stream));

  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaInDegree, 0,
                                  deviceGraph.numNodes * sizeof(int),
                                  stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaOutDegree, 0,
                                  deviceGraph.numNodes * sizeof(int),
                                  stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaComponentSizes, 0,
                                  numComponents * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaBadComponents, 0,
                                  numComponents * sizeof(int), stream));

  const dim3 edgeGrid((deviceGraph.numEdges + kBlockSize - 1) / kBlockSize);
  const dim3 nodeGrid((deviceGraph.numNodes + kBlockSize - 1) / kBlockSize);
  computeComponentDegreeKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numEdges, deviceGraph.src, deviceGraph.dst,
      deviceGraph.valid, cudaInDegree, cudaOutDegree);
  ACTS_CUDA_CHECK(cudaGetLastError());
  countActiveComponentNodesKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numNodes, deviceGraph.activeNodes, cudaLabels,
      cudaComponentSizes);
  ACTS_CUDA_CHECK(cudaGetLastError());
  markBadComponentsKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numNodes, deviceGraph.activeNodes, cudaLabels, cudaInDegree,
      cudaOutDegree, cudaBadComponents);
  ACTS_CUDA_CHECK(cudaGetLastError());
  buildInitialComponentMasksKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      deviceGraph.numNodes, deviceGraph.activeNodes, cudaLabels,
      cudaComponentSizes, cudaBadComponents, minCandidateSize,
      *cudaSimpleNodeMask, *cudaComplexNodeMask);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(cudaInDegree, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaOutDegree, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaComponentSizes, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaBadComponents, stream));
}

DeviceCsrGraph buildCsrGraphCuda(const DeviceCompactEdges &edges,
                                 const std::string &pathMetric,
                                 cudaStream_t stream) {
  DeviceCsrGraph graph;
  graph.numNodes = edges.numNodes;
  graph.numEdges = edges.numEdges;
  if (edges.numEdges == 0) {
    return graph;
  }

  int *sortedSrc{};
  int *sortedDst{};
  float *sortedScore{};
  int *sortedIncomingDst{};
  int *sortedIncomingSrc{};
  int *rowCounts{};
  int *incomingCounts{};
  int *cursor{};
  int *incomingCursor{};

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&sortedSrc, edges.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&sortedDst, edges.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&sortedScore, edges.numEdges * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&sortedIncomingDst,
                                  edges.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&sortedIncomingSrc,
                                  edges.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&rowCounts, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&incomingCounts, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cursor, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&incomingCursor, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&graph.rowPtr, (edges.numNodes + 1) * sizeof(int),
                      stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&graph.colIdx, edges.numEdges * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&graph.edgeWeight,
                                  edges.numEdges * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&graph.incomingRowPtr,
                                  (edges.numNodes + 1) * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&graph.incomingColIdx,
                                  edges.numEdges * sizeof(int), stream));

  ACTS_CUDA_CHECK(
      cudaMemsetAsync(rowCounts, 0, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(incomingCounts, 0,
                                  edges.numNodes * sizeof(int), stream));
  const dim3 edgeGrid((edges.numEdges + kBlockSize - 1) / kBlockSize);
  fillCountsKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      edges.numEdges, edges.src, rowCounts);
  fillCountsKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      edges.numEdges, edges.dst, incomingCounts);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaMemsetAsync(graph.rowPtr, 0, sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(graph.incomingRowPtr, 0, sizeof(int),
                                  stream));
  thrust::inclusive_scan(thrust::device.on(stream), rowCounts,
                         rowCounts + edges.numNodes, graph.rowPtr + 1);
  thrust::inclusive_scan(thrust::device.on(stream), incomingCounts,
                         incomingCounts + edges.numNodes,
                         graph.incomingRowPtr + 1);

  ACTS_CUDA_CHECK(cudaMemcpyAsync(sortedSrc, edges.src,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(sortedDst, edges.dst,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(sortedScore, edges.score,
                                  edges.numEdges * sizeof(float),
                                  cudaMemcpyDeviceToDevice, stream));
  thrust::sort_by_key(thrust::device.on(stream), sortedSrc,
                      sortedSrc + edges.numEdges,
                      thrust::make_zip_iterator(
                          thrust::make_tuple(sortedDst, sortedScore)));

  ACTS_CUDA_CHECK(cudaMemcpyAsync(sortedIncomingDst, edges.dst,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(sortedIncomingSrc, edges.src,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToDevice, stream));
  thrust::sort_by_key(thrust::device.on(stream), sortedIncomingDst,
                      sortedIncomingDst + edges.numEdges, sortedIncomingSrc);

  ACTS_CUDA_CHECK(
      cudaMemsetAsync(cursor, 0, edges.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(incomingCursor, 0,
                                  edges.numNodes * sizeof(int), stream));
  scatterSortedOutgoingKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      edges.numEdges, sortedSrc, sortedDst, sortedScore, graph.rowPtr, cursor,
      graph.colIdx, graph.edgeWeight, pathMetric == "length");
  scatterSortedIncomingKernel<<<edgeGrid, kBlockSize, 0, stream>>>(
      edges.numEdges, sortedIncomingDst, sortedIncomingSrc,
      graph.incomingRowPtr, incomingCursor, graph.incomingColIdx);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(sortedSrc, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(sortedDst, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(sortedScore, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(sortedIncomingDst, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(sortedIncomingSrc, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(rowCounts, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(incomingCounts, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cursor, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(incomingCursor, stream));
  return graph;
}

DpCudaState runDpOnCsrCuda(const DeviceCsrGraph &graph,
                           const unsigned char *cudaActiveNodes,
                           cudaStream_t stream) {
  DpCudaState state;
  state.numNodes = graph.numNodes;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&state.bestScore, graph.numNodes * sizeof(float),
                      stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&state.bestChild, graph.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&state.sourceMask,
                                  graph.numNodes * sizeof(unsigned char),
                                  stream));
  if (graph.numEdges == 0) {
    return state;
  }

  int *cudaInDegree{};
  int *cudaRemainingOutDegree{};
  int *cudaFrontier{};
  int *cudaNextFrontier{};
  int *cudaFrontierSize{};
  int *cudaNextFrontierSize{};
  float *cudaFrontierBestScore{};
  int *cudaFrontierBestChild{};

  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaInDegree,
                                  graph.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaRemainingOutDegree,
                                  graph.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaFrontier, graph.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaNextFrontier, graph.numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaFrontierSize, sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaNextFrontierSize, sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaFrontierBestScore,
                                  graph.numNodes * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaFrontierBestChild,
                                  graph.numNodes * sizeof(int), stream));

  const dim3 nodeGrid((graph.numNodes + kBlockSize - 1) / kBlockSize);
  initializeDpKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      state.bestScore, state.bestChild, state.sourceMask, cudaInDegree,
      cudaRemainingOutDegree, graph.numNodes);
  ACTS_CUDA_CHECK(cudaGetLastError());
  computeActiveDegreeKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      graph.rowPtr, graph.colIdx, cudaActiveNodes, cudaInDegree,
      cudaRemainingOutDegree, graph.numNodes);
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaFrontierSize, 0, sizeof(int), stream));
  initializeFrontierKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      cudaActiveNodes, cudaRemainingOutDegree, cudaFrontier, cudaFrontierSize,
      graph.numNodes);
  ACTS_CUDA_CHECK(cudaGetLastError());

  int currentFrontierSize = 0;
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&currentFrontierSize, cudaFrontierSize,
                                  sizeof(int), cudaMemcpyDeviceToHost,
                                  stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  while (currentFrontierSize > 0) {
    const dim3 frontierGrid((currentFrontierSize + kBlockSize - 1) /
                            kBlockSize);
    processFrontierKernel<<<frontierGrid, kBlockSize, 0, stream>>>(
        cudaFrontier, currentFrontierSize, graph.rowPtr, graph.colIdx,
        graph.edgeWeight, cudaActiveNodes, state.bestScore,
        cudaFrontierBestScore, cudaFrontierBestChild);
    ACTS_CUDA_CHECK(cudaGetLastError());
    finalizeFrontierKernel<<<frontierGrid, kBlockSize, 0, stream>>>(
        cudaFrontier, currentFrontierSize, cudaFrontierBestScore,
        cudaFrontierBestChild, state.bestScore, state.bestChild);
    ACTS_CUDA_CHECK(cudaGetLastError());
    ACTS_CUDA_CHECK(
        cudaMemsetAsync(cudaNextFrontierSize, 0, sizeof(int), stream));
    enqueueParentFrontierKernel<<<frontierGrid, kBlockSize, 0, stream>>>(
        cudaFrontier, currentFrontierSize, graph.incomingRowPtr,
        graph.incomingColIdx, cudaActiveNodes, cudaRemainingOutDegree,
        cudaNextFrontier, cudaNextFrontierSize);
    ACTS_CUDA_CHECK(cudaGetLastError());
    std::swap(cudaFrontier, cudaNextFrontier);
    std::swap(cudaFrontierSize, cudaNextFrontierSize);
    ACTS_CUDA_CHECK(cudaMemcpyAsync(&currentFrontierSize, cudaFrontierSize,
                                    sizeof(int), cudaMemcpyDeviceToHost,
                                    stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  buildSourceMaskKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      cudaActiveNodes, cudaInDegree, state.bestScore, state.sourceMask,
      graph.numNodes);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(cudaInDegree, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaRemainingOutDegree, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaFrontier, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaNextFrontier, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaFrontierSize, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaNextFrontierSize, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaFrontierBestScore, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaFrontierBestChild, stream));
  return state;
}

std::pair<std::vector<int>, std::vector<int>> selectAndTracePathsCuda(
    const DpCudaState &dpState, const int *cudaComponentLabels,
    int numComponents, float minRootScore, unsigned char *cudaActiveNodes,
    cudaStream_t stream) {
  if (numComponents <= 0) {
    return {};
  }

  const std::size_t numNodes = dpState.numNodes;
  int *cudaSelectedRoots{};
  float *cudaComponentBestScore{};
  int *cudaSelectedTrackLabels{};
  int *cudaSelectedNodes{};
  int *cudaSelectedCount{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaSelectedRoots,
                                  numComponents * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaComponentBestScore,
                                  numComponents * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaSelectedTrackLabels,
                                  numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaSelectedNodes, numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaSelectedCount, sizeof(int), stream));

  const dim3 componentGrid((numComponents + kBlockSize - 1) / kBlockSize);
  initFloatKernel<<<componentGrid, kBlockSize, 0, stream>>>(
      numComponents, -std::numeric_limits<float>::infinity(),
      cudaComponentBestScore);
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(cudaSelectedRoots, 0xff, numComponents * sizeof(int),
                      stream));
  const dim3 nodeGrid((numNodes + kBlockSize - 1) / kBlockSize);
  selectComponentScoresKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      numNodes, cudaComponentLabels, dpState.sourceMask, dpState.bestScore,
      cudaComponentBestScore);
  ACTS_CUDA_CHECK(cudaGetLastError());
  selectRootsKernel<<<nodeGrid, kBlockSize, 0, stream>>>(
      numNodes, cudaComponentLabels, dpState.sourceMask, dpState.bestScore,
      cudaComponentBestScore, minRootScore, cudaSelectedRoots);
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaSelectedCount, 0, sizeof(int), stream));

  const dim3 pathGrid((numComponents + kBlockSize - 1) / kBlockSize);
  traceSelectedPathsKernel<<<pathGrid, kBlockSize, 0, stream>>>(
      dpState.bestChild, cudaSelectedRoots, cudaSelectedTrackLabels,
      cudaSelectedNodes, cudaSelectedCount, numNodes, numComponents);
  ACTS_CUDA_CHECK(cudaGetLastError());

  int selectedCount = 0;
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&selectedCount, cudaSelectedCount,
                                  sizeof(int), cudaMemcpyDeviceToHost,
                                  stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  if (selectedCount > 0 && cudaActiveNodes != nullptr) {
    const dim3 selectedGrid((selectedCount + kBlockSize - 1) / kBlockSize);
    deactivateSelectedNodesKernel<<<selectedGrid, kBlockSize, 0, stream>>>(
        selectedCount, cudaSelectedNodes, cudaActiveNodes);
    ACTS_CUDA_CHECK(cudaGetLastError());
  }

  std::vector<int> labels(selectedCount);
  std::vector<int> nodes(selectedCount);
  if (selectedCount > 0) {
    ACTS_CUDA_CHECK(cudaMemcpyAsync(labels.data(), cudaSelectedTrackLabels,
                                    selectedCount * sizeof(int),
                                    cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaMemcpyAsync(nodes.data(), cudaSelectedNodes,
                                    selectedCount * sizeof(int),
                                    cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSelectedRoots, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaComponentBestScore, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSelectedTrackLabels, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSelectedNodes, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSelectedCount, stream));

  return {std::move(labels), std::move(nodes)};
}

void appendSmallResidualDWalkTracks(
    const DeviceCompactEdges &edges,
    const std::vector<int> &complexSpacePointIds,
    std::size_t minCandidateSize, const std::string &pathMetric,
    cudaStream_t stream, std::vector<std::vector<int>> &trackCandidates) {
  if (edges.numEdges == 0) {
    return;
  }

  std::vector<int> src(edges.numEdges);
  std::vector<int> dst(edges.numEdges);
  std::vector<float> score(edges.numEdges);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(src.data(), edges.src,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(dst.data(), edges.dst,
                                  edges.numEdges * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(score.data(), edges.score,
                                  edges.numEdges * sizeof(float),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  std::vector<unsigned char> activeEdge(edges.numEdges, 1);
  const float rootScoreCut = minRootScore(pathMetric);
  const bool useLengthMetric = pathMetric == "length";

  while (true) {
    std::vector<std::size_t> currentEdges;
    currentEdges.reserve(edges.numEdges);
    std::vector<unsigned char> activeNode(edges.numNodes, 0);
    for (std::size_t edge = 0; edge < edges.numEdges; ++edge) {
      if (activeEdge.at(edge) == 0) {
        continue;
      }
      currentEdges.push_back(edge);
      activeNode.at(src.at(edge)) = 1;
      activeNode.at(dst.at(edge)) = 1;
    }
    if (currentEdges.empty()) {
      break;
    }

    DisjointSet components(edges.numNodes);
    std::vector<int> inDegree(edges.numNodes, 0);
    std::vector<std::vector<std::pair<int, float>>> adjacency(edges.numNodes);
    for (std::size_t edge : currentEdges) {
      components.unite(src.at(edge), dst.at(edge));
      ++inDegree.at(dst.at(edge));
      adjacency.at(src.at(edge))
          .push_back({dst.at(edge), useLengthMetric ? 1.0F : score.at(edge)});
    }

    std::unordered_map<int, int> componentMap;
    std::vector<int> componentLabel(edges.numNodes, -1);
    int numComponents = 0;
    for (std::size_t node = 0; node < edges.numNodes; ++node) {
      if (activeNode.at(node) == 0) {
        continue;
      }
      auto [it, inserted] =
          componentMap.emplace(components.find(static_cast<int>(node)),
                               numComponents);
      if (inserted) {
        ++numComponents;
      }
      componentLabel.at(node) = it->second;
    }

    std::vector<float> bestScore(edges.numNodes, 0.0F);
    std::vector<int> bestChild(edges.numNodes, -1);
    std::vector<unsigned char> visitState(edges.numNodes, 0);
    std::function<float(int)> solve = [&](int node) -> float {
      if (visitState.at(node) == 2) {
        return bestScore.at(node);
      }
      if (visitState.at(node) == 1) {
        return 0.0F;
      }
      visitState.at(node) = 1;
      float bestValue = 0.0F;
      int bestNext = -1;
      for (auto [child, weight] : adjacency.at(node)) {
        if (activeNode.at(child) == 0) {
          continue;
        }
        float candidate = weight + solve(child);
        if (candidate > bestValue) {
          bestValue = candidate;
          bestNext = child;
        }
      }
      bestScore.at(node) = bestValue;
      bestChild.at(node) = bestNext;
      visitState.at(node) = 2;
      return bestValue;
    };

    std::vector<float> componentBest(numComponents,
                                     -std::numeric_limits<float>::infinity());
    std::vector<int> componentRoot(numComponents, -1);
    for (std::size_t node = 0; node < edges.numNodes; ++node) {
      if (activeNode.at(node) == 0 || inDegree.at(node) != 0) {
        continue;
      }
      float value = solve(static_cast<int>(node));
      int component = componentLabel.at(node);
      if (component >= 0 && value > rootScoreCut &&
          value > componentBest.at(component)) {
        componentBest.at(component) = value;
        componentRoot.at(component) = static_cast<int>(node);
      }
    }

    std::vector<unsigned char> selectedNode(edges.numNodes, 0);
    bool selectedAny = false;
    for (int root : componentRoot) {
      if (root < 0) {
        continue;
      }
      std::vector<int> path;
      int node = root;
      std::size_t guard = 0;
      while (node >= 0 && guard < edges.numNodes &&
             selectedNode.at(node) == 0) {
        path.push_back(node);
        selectedNode.at(node) = 1;
        selectedAny = true;
        node = bestChild.at(node);
        ++guard;
      }
      if (path.size() < minCandidateSize) {
        continue;
      }
      std::vector<int> track;
      track.reserve(path.size());
      for (int pathNode : path) {
        track.push_back(complexSpacePointIds.at(pathNode));
      }
      trackCandidates.push_back(std::move(track));
    }

    if (!selectedAny) {
      break;
    }
    for (std::size_t edge = 0; edge < edges.numEdges; ++edge) {
      if (activeEdge.at(edge) != 0 &&
          (selectedNode.at(src.at(edge)) != 0 ||
           selectedNode.at(dst.at(edge)) != 0)) {
        activeEdge.at(edge) = 0;
      }
    }
  }
}

}  // namespace

namespace ActsPlugins {

std::vector<std::vector<int>> DWalkTrackBuilding::operator()(
    PipelineTensors tensors, std::vector<int> &spacePointIds,
    const ExecutionContext &execContext) {
  ACTS_VERBOSE("Start CUDA D-WALK track building");

  if (!tensors.edgeScores.has_value()) {
    throw std::runtime_error("DWalkTrackBuilding expects edge scores");
  }
  if (!(tensors.edgeIndex.device().isCuda() &&
        tensors.edgeScores->device().isCuda() &&
        tensors.nodeFeatures.device().isCuda())) {
    throw std::runtime_error(
        "DWalkTrackBuilding expects tensors to be on CUDA");
  }
  if (m_cfg.pathMetric != "score_weighted_length" &&
      m_cfg.pathMetric != "length") {
    throw std::invalid_argument(
        "DWalkTrackBuilding pathMetric must be 'score_weighted_length' or "
        "'length'");
  }

  assert(tensors.edgeIndex.shape().at(0) == 2);
  assert(tensors.edgeIndex.shape().at(1) == tensors.edgeScores->shape().at(0));

  const auto numNodes = tensors.nodeFeatures.shape().at(0);
  const auto numFeatures = tensors.nodeFeatures.shape().at(1);
  const auto numEdges = tensors.edgeIndex.shape().at(1);
  if (m_cfg.radialFeatureIndex >= numFeatures) {
    throw std::out_of_range("DWalkTrackBuilding radialFeatureIndex is invalid");
  }
  if (numNodes > spacePointIds.size()) {
    throw std::runtime_error(
        "DWalkTrackBuilding received more graph nodes than space point IDs");
  }
  if (numEdges == 0 || numNodes < m_cfg.minCandidateSize) {
    return {};
  }

  auto stream = execContext.stream.value();
  std::vector<int> initialLabels;
  DeviceOrientedEdges deviceGraph;
  int *cudaInitialLabels{};
  int initialNumComponents = 0;
  auto edges = createOrientedEdgesCuda(tensors.edgeIndex, *tensors.edgeScores,
                                       tensors.nodeFeatures,
                                       m_cfg.radialFeatureIndex, stream,
                                       deviceGraph, &cudaInitialLabels,
                                       &initialNumComponents, &initialLabels);
  if (edges.empty()) {
    freeDeviceOrientedEdges(deviceGraph, stream);
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaInitialLabels, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    return {};
  }

  unsigned char *cudaSimpleNodeMask{};
  unsigned char *cudaComplexNodeMask{};
  classifyInitialComponentsCuda(
      deviceGraph, cudaInitialLabels, initialNumComponents,
      static_cast<int>(m_cfg.minCandidateSize), &cudaSimpleNodeMask,
      &cudaComplexNodeMask, stream);

  std::vector<unsigned char> simpleNodeMask(numNodes, 0);
  std::vector<unsigned char> complexNodeMask(numNodes, 0);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(simpleNodeMask.data(), cudaSimpleNodeMask,
                                  numNodes * sizeof(unsigned char),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(complexNodeMask.data(), cudaComplexNodeMask,
                                  numNodes * sizeof(unsigned char),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  std::vector<int> inDegree(numNodes, 0);
  std::vector<int> outDegree(numNodes, 0);
  std::vector<int> nextNode(numNodes, -1);
  for (const auto &edge : edges) {
    ++outDegree.at(edge.src);
    ++inDegree.at(edge.dst);
    nextNode.at(edge.src) = edge.dst;
  }

  std::unordered_map<int, std::vector<int>> componentNodes;
  for (std::size_t node = 0; node < numNodes; ++node) {
    if (simpleNodeMask.at(node) != 0) {
      const int label = initialLabels.at(node);
      componentNodes[label].push_back(static_cast<int>(node));
    }
  }
  std::vector<std::vector<int>> trackCandidates;
  for (const auto &[component, nodes] : componentNodes) {
    auto orderedNodes = orderedSimpleComponentNodes(
        nodes, component, initialLabels, inDegree, nextNode);
    std::vector<int> track;
    track.reserve(orderedNodes.size());
    for (int node : orderedNodes) {
      track.push_back(spacePointIds.at(node));
    }
    trackCandidates.push_back(std::move(track));
  }

  std::vector<int> originalToComplex(numNodes, -1);
  std::vector<int> complexSpacePointIds;
  for (std::size_t node = 0; node < numNodes; ++node) {
    if (complexNodeMask.at(node) != 0) {
      originalToComplex.at(node) = static_cast<int>(complexSpacePointIds.size());
      complexSpacePointIds.push_back(spacePointIds.at(node));
    }
  }
  const std::size_t numComplexNodes = complexSpacePointIds.size();
  ACTS_DEBUG("CUDA D-WALK complex nodes: " << numComplexNodes);

  if (numComplexNodes == 0) {
    freeDeviceOrientedEdges(deviceGraph, stream);
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaInitialLabels, stream));
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaSimpleNodeMask, stream));
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaComplexNodeMask, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("CUDA D-WALK found " << trackCandidates.size()
                                    << " track candidates");
    return trackCandidates;
  }

  int *cudaOriginalToComplex{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaOriginalToComplex,
                                  numNodes * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaOriginalToComplex,
                                  originalToComplex.data(),
                                  numNodes * sizeof(int),
                                  cudaMemcpyHostToDevice, stream));

  auto complexEdges = maxAddCompactDeviceEdgesCuda(
      deviceGraph, cudaComplexNodeMask, m_cfg.thMin, m_cfg.thAdd,
      stream);
  freeDeviceOrientedEdges(deviceGraph, stream);
  if (complexEdges.numEdges != 0) {
    const dim3 remapGrid((complexEdges.numEdges + kBlockSize - 1) / kBlockSize);
    remapEdgesKernel<<<remapGrid, kBlockSize, 0, stream>>>(
        complexEdges.numEdges, complexEdges.src, complexEdges.dst,
        cudaOriginalToComplex);
    ACTS_CUDA_CHECK(cudaGetLastError());
    complexEdges.numNodes = numComplexNodes;
  }
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaInitialLabels, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSimpleNodeMask, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaComplexNodeMask, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaOriginalToComplex, stream));
  ACTS_DEBUG("CUDA D-WALK initial complex edges: " << complexEdges.numEdges);

  unsigned char *cudaComplexActiveNodes{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaComplexActiveNodes,
                                  numComplexNodes * sizeof(unsigned char),
                                  stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaComplexActiveNodes, 1,
                                  numComplexNodes * sizeof(unsigned char),
                                  stream));

  int iteration = 0;
  while (complexEdges.numEdges != 0) {
    if (complexEdges.numEdges <= kSmallResidualEdgeThreshold) {
      appendSmallResidualDWalkTracks(complexEdges, complexSpacePointIds,
                                     m_cfg.minCandidateSize, m_cfg.pathMetric,
                                     stream, trackCandidates);
      break;
    }

    ++iteration;
    int *cudaCurrentLabels{};
    ACTS_CUDA_CHECK(
        cudaMallocAsync(&cudaCurrentLabels,
                        numComplexNodes * sizeof(int), stream));
    int currentNumComponents = ActsPlugins::detail::connectedComponentsCuda(
        complexEdges.numEdges, complexEdges.src, complexEdges.dst,
        numComplexNodes, cudaCurrentLabels, stream, false);
    ACTS_DEBUG("CUDA D-WALK iteration " << iteration
               << ": components=" << currentNumComponents
               << ", edges=" << complexEdges.numEdges);
    if (currentNumComponents == 0) {
      ACTS_CUDA_CHECK(cudaFreeAsync(cudaCurrentLabels, stream));
      break;
    }

    auto csrGraph = buildCsrGraphCuda(complexEdges, m_cfg.pathMetric, stream);
    auto dpState = runDpOnCsrCuda(csrGraph, cudaComplexActiveNodes, stream);
    freeDeviceCsrGraph(csrGraph, stream);

    auto [selectedTrackLabels, selectedNodes] =
        selectAndTracePathsCuda(dpState, cudaCurrentLabels,
                                currentNumComponents,
                                minRootScore(m_cfg.pathMetric),
                                cudaComplexActiveNodes, stream);
    freeDpCudaState(dpState, stream);
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaCurrentLabels, stream));
    if (selectedNodes.empty()) {
      break;
    }

    std::vector<std::vector<int>> paths(currentNumComponents);
    for (std::size_t selected = 0; selected < selectedNodes.size(); ++selected) {
      int label = selectedTrackLabels.at(selected);
      int node = selectedNodes.at(selected);
      if (label < 0 || label >= static_cast<int>(paths.size()) || node < 0 ||
          node >= static_cast<int>(numComplexNodes)) {
        continue;
      }
      paths.at(label).push_back(node);
    }

    for (const auto &path : paths) {
      if (path.size() < m_cfg.minCandidateSize) {
        continue;
      }
      std::vector<int> track;
      track.reserve(path.size());
      for (int node : path) {
        track.push_back(complexSpacePointIds.at(node));
      }
      trackCandidates.push_back(std::move(track));
    }

    auto updatedEdges = compactActiveEdgesCuda(complexEdges,
                                               cudaComplexActiveNodes, stream);
    freeDeviceCompactEdges(complexEdges, stream);
    complexEdges = updatedEdges;
  }

  freeDeviceCompactEdges(complexEdges, stream);
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaComplexActiveNodes, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ACTS_DEBUG("CUDA D-WALK found " << trackCandidates.size()
                                  << " track candidates");
  return trackCandidates;
}

}  // namespace ActsPlugins
