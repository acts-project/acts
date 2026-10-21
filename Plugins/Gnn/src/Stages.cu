// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/Stages.hpp"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

namespace {

/// Scatter true into mask[usedNodes[i]] for each i in [0, nUsed)
__global__ void buildNodeMaskKernel(const std::int64_t *usedNodes,
                                    std::size_t nUsed, bool *mask) {
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nUsed) {
    mask[usedNodes[i]] = true;
  }
}

/// Scatter new index newIdx into oldToNew[usedNodes[newIdx]] for each newIdx
__global__ void buildReverseMapKernel(const std::int64_t *usedNodes,
                                      std::size_t nUsed,
                                      std::int64_t *oldToNew) {
  const std::size_t newIdx = blockIdx.x * blockDim.x + threadIdx.x;
  if (newIdx < nUsed) {
    oldToNew[usedNodes[newIdx]] = static_cast<std::int64_t>(newIdx);
  }
}

/// Apply old→new remapping to all 2*nEdges edge endpoint values in-place
__global__ void remapEdgesKernel(std::size_t nTotal,
                                 const std::int64_t *oldToNew,
                                 std::int64_t *edgeData) {
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nTotal) {
    edgeData[i] = oldToNew[edgeData[i]];
  }
}

}  // namespace

namespace ActsPlugins::detail {

PipelineTensors cudaRemoveUnusedNodes(PipelineTensors &&tensors,
                                      std::vector<int> &spacePointIds,
                                      const ExecutionContext &execCtx) {
  const auto stream = execCtx.stream.value();
  const auto nNodes = tensors.nodeFeatures.shape()[0];
  const auto nEdges = tensors.edgeIndex.shape()[1];

  // Copy edgeIndex into a scratch buffer — thrust needs a mutable working
  // copy and the edgeIndex tensor is remapped in-place later.
  auto tmp = Tensor<std::int64_t>::Create({1, 2 * nEdges}, execCtx);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(tmp.data(), tensors.edgeIndex.data(),
                                  tmp.nbytes(), cudaMemcpyDeviceToDevice,
                                  stream));

  // Sort + unique → sorted unique used-node indices in tmp[0..nUsed)
  thrust::sort(thrust::device.on(stream), tmp.data(), tmp.data() + 2 * nEdges);
  auto *uniqEnd = thrust::unique(thrust::device.on(stream), tmp.data(),
                                 tmp.data() + 2 * nEdges);
  // nUsed must be read on host — sync the stream just for this scalar
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  const std::size_t nUsed = static_cast<std::size_t>(uniqEnd - tmp.data());

  // Boolean node mask [nNodes, 1]: true at each surviving node index
  auto mask = Tensor<bool>::Create({nNodes, 1}, execCtx);
  ACTS_CUDA_CHECK(cudaMemsetAsync(mask.data(), 0, mask.nbytes(), stream));
  const dim3 blockDim = 1024;
  const dim3 gridUsed = (nUsed + blockDim.x - 1) / blockDim.x;
  buildNodeMaskKernel<<<gridUsed, blockDim, 0, stream>>>(tmp.data(), nUsed,
                                                         mask.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Reverse map [nNodes, 1]: oldToNew[old] = new index after compaction
  auto oldToNew = Tensor<std::int64_t>::Create({nNodes, 1}, execCtx);
  buildReverseMapKernel<<<gridUsed, blockDim, 0, stream>>>(tmp.data(), nUsed,
                                                           oldToNew.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Shrink nodeFeatures to surviving rows using the mask
  auto newNodeFeatures = selectRows(tensors.nodeFeatures, mask, execCtx);

  // Remap edge endpoint indices in-place using the old→new map
  const dim3 gridEdges = (2 * nEdges + blockDim.x - 1) / blockDim.x;
  remapEdgesKernel<<<gridEdges, blockDim, 0, stream>>>(
      2 * nEdges, oldToNew.data(), tensors.edgeIndex.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Copy surviving node indices to host to update spacePointIds
  std::vector<std::int64_t> hostUsedNodes(nUsed);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(hostUsedNodes.data(), tmp.data(),
                                  nUsed * sizeof(std::int64_t),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  std::vector<int> remapped;
  remapped.reserve(nUsed);
  for (const auto oldIdx : hostUsedNodes) {
    remapped.push_back(spacePointIds[static_cast<std::size_t>(oldIdx)]);
  }
  spacePointIds = std::move(remapped);

  return {std::move(newNodeFeatures), std::move(tensors.edgeIndex),
          std::move(tensors.edgeFeatures), std::move(tensors.edgeScores)};
}

}  // namespace ActsPlugins::detail
