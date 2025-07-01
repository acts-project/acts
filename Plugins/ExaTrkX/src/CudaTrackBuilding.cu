// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/CudaTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/detail/ConnectedComponents.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"
#include "Acts/Plugins/ExaTrkX/detail/JunctionRemoval.hpp"
#include "Acts/Utilities/Zip.hpp"

namespace Acts {

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

  // TODO not sure why there is an issue that is not detected in the unit tests
  numberLabels += 1;

  std::vector<int> trackLabels(numSpacepoints);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(trackLabels.data(), cudaLabels,
                                  numSpacepoints * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));

  // Free Memory
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaLabels, stream));
  if (m_cfg.doJunctionRemoval) {
    ACTS_CUDA_CHECK(cudaFreeAsync(cudaSrcPtr, stream));
  }

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("Found " << numberLabels << " track candidates");

  std::vector<std::vector<int>> trackCandidates(numberLabels);

  for (const auto [label, id] : Acts::zip(trackLabels, spacepointIDs)) {
    trackCandidates[label].push_back(id);
  }

  return trackCandidates;
}

}  // namespace Acts
