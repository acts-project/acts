// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Zip.hpp"
#include "ActsPlugins/Gnn/CudaTrackBuilding.hpp"
#include "ActsPlugins/Gnn/detail/ConnectedComponents.cuh"
#include "ActsPlugins/Gnn/detail/CudaUtils.cuh"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"
#include "ActsPlugins/Gnn/detail/JunctionRemoval.hpp"

using namespace Acts;

namespace ActsPlugins {

std::vector<std::vector<int>> CudaTrackBuilding::operator()(
    PipelineTensors tensors, std::vector<int>& spacepointIDs,
    const ExecutionContext& execContext) {
  ACTS_VERBOSE("Start CUDA track building");
  if (!(tensors.edgeIndex.device().isCuda() &&
        tensors.edgeScores.value().device().isCuda())) {
    throw std::runtime_error(
        "CudaTrackBuilding expects tensors to be on CUDA!");
  }

  const auto numSpacepoints = spacepointIDs.size();
  auto numEdges = static_cast<std::size_t>(tensors.edgeIndex.shape().at(1));

  if (numEdges == 0) {
    ACTS_DEBUG("No edges remained after edge classification");
    return {};
  }

  auto stream = execContext.stream.value();

  auto cudaSrcPtr = tensors.edgeIndex.data();
  auto cudaTgtPtr = tensors.edgeIndex.data() + numEdges;

  auto ms = [](auto t0, auto t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
        .count();
  };

  if (m_cfg.doJunctionRemoval) {
    assert(tensors.edgeScores->shape().at(0) ==
           tensors.edgeIndex.shape().at(1));
    auto cudaScorePtr = tensors.edgeScores->data();

    ACTS_DEBUG("Do junction removal...");
    auto t0 = std::chrono::high_resolution_clock::now();
    auto [cudaSrcPtrJr, numEdgesOut] = detail::junctionRemovalCuda(
        numEdges, numSpacepoints, cudaScorePtr, cudaSrcPtr, cudaTgtPtr, stream);
    auto t1 = std::chrono::high_resolution_clock::now();
    cudaSrcPtr = cudaSrcPtrJr;
    cudaTgtPtr = cudaSrcPtrJr + numEdgesOut;

    if (numEdgesOut == 0) {
      ACTS_WARNING(
          "No edges remained after junction removal, this should not happen!");
      ACTS_CUDA_CHECK(cudaFreeAsync(cudaSrcPtrJr, stream));
      ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
      return {};
    }

    ACTS_DEBUG("Removed " << numEdges - numEdgesOut
                          << " edges in junction removal");
    ACTS_DEBUG("Junction removal took " << ms(t0, t1) << " ms");
    numEdges = numEdgesOut;
  }

  int* cudaLabels{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaLabels, numSpacepoints * sizeof(int), stream));

  auto t0 = std::chrono::high_resolution_clock::now();
  std::size_t numberLabels = detail::connectedComponentsCuda(
      numEdges, cudaSrcPtr, cudaTgtPtr, numSpacepoints, cudaLabels, stream,
      m_cfg.useOneBlockImplementation);
  auto t1 = std::chrono::high_resolution_clock::now();
  ACTS_DEBUG("Connected components took " << ms(t0, t1) << " ms");
  ACTS_VERBOSE("Found " << numberLabels << " track candidates");

  // Postprocess labels
  int* cudaSpacepointIDs{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaSpacepointIDs,
                                  spacepointIDs.size() * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaSpacepointIDs, spacepointIDs.data(),
                                  spacepointIDs.size() * sizeof(int),
                                  cudaMemcpyHostToDevice, stream));

  // Allocate space for the bounds
  int* cudaBounds{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaBounds, (numberLabels + 1) * sizeof(int), stream));

  // Compute the bounds of the track candidates
  detail::findTrackCandidateBounds(cudaLabels, cudaSpacepointIDs, cudaBounds,
                                   numSpacepoints, numberLabels, stream);

  // Copy the bounds to the host
  std::vector<int> bounds(numberLabels + 1);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(bounds.data(), cudaBounds,
                                  (numberLabels + 1) * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));

  // Copy the sorted spacepoint IDs to the host
  ACTS_CUDA_CHECK(cudaMemcpyAsync(spacepointIDs.data(), cudaSpacepointIDs,
                                  spacepointIDs.size() * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));

  // Free Memory
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaLabels, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaSpacepointIDs, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaBounds, stream));
  if (m_cfg.doJunctionRemoval) {
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaSrcPtr, stream));
  }

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_DEBUG("Bounds size: " << bounds.size());
  ACTS_DEBUG("Bounds: " << bounds.at(0) << ", " << bounds.at(1) << ", "
                        << bounds.at(2) << ", ..., "
                        << bounds.at(numberLabels - 2) << ", "
                        << bounds.at(numberLabels - 1) << ", "
                        << bounds.at(numberLabels));

  std::vector<std::vector<int>> trackCandidates;
  trackCandidates.reserve(numberLabels);
  for (std::size_t label = 0ul; label < numberLabels; ++label) {
    int start = bounds.at(label);
    int end = bounds.at(label + 1);

    assert(start >= 0);
    assert(end <= static_cast<int>(numSpacepoints));
    assert(start <= end);

    if (end - start < m_cfg.minCandidateSize) {
      continue;
    }

    trackCandidates.emplace_back(spacepointIDs.begin() + start,
                                 spacepointIDs.begin() + end);
  }

  return trackCandidates;
}

}  // namespace ActsPlugins
