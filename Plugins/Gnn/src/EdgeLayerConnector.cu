// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/EdgeLayerConnector.hpp"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <MMG/CUDA_edge_layer_connector>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>

using namespace Acts;

namespace ActsPlugins {

std::vector<std::vector<int>> EdgeLayerConnector::operator()(
    PipelineTensors tensors, std::vector<int>& spacepointIDs,
    const ExecutionContext& execContext) {
  ACTS_DEBUG("Get event data");

  const auto numEdges = static_cast<int>(tensors.edgeIndex.shape().at(1));
  const auto numSpacepoints = static_cast<int>(spacepointIDs.size());

  auto stream = execContext.stream.value();

  // Convert std::int64_t edge indices to int using Tensor for memory management
  auto srcInt64Ptr = tensors.edgeIndex.data();
  auto tgtInt64Ptr = tensors.edgeIndex.data() + numEdges;

  auto edgeSrc =
      Tensor<int>::Create({1, static_cast<std::size_t>(numEdges)}, execContext);
  auto edgeTgt =
      Tensor<int>::Create({1, static_cast<std::size_t>(numEdges)}, execContext);

  thrust::copy(thrust::cuda::par.on(stream), srcInt64Ptr,
               srcInt64Ptr + numEdges, edgeSrc.data());
  thrust::copy(thrust::cuda::par.on(stream), tgtInt64Ptr,
               tgtInt64Ptr + numEdges, edgeTgt.data());

  // Copy spacepoint IDs to GPU
  auto spacepointIDsTensor =
      Tensor<int>::Create({1, spacepointIDs.size()}, execContext);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(
      spacepointIDsTensor.data(), spacepointIDs.data(),
      spacepointIDs.size() * sizeof(int), cudaMemcpyHostToDevice, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ACTS_DEBUG("Setup graph...");
  CUDA_graph<float> graph(spacepointIDsTensor.data(), numSpacepoints,
                          edgeSrc.data(), edgeTgt.data(),
                          tensors.edgeScores->data(), numEdges);
  ACTS_CUDA_CHECK(cudaDeviceSynchronize());
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_DEBUG("Setup EdgeLayerConnector...");
  CUDA_edge_layer_connector<float> connector(
      &graph, m_cfg.weightsCut, static_cast<int>(m_cfg.minHits),
      static_cast<int>(m_cfg.blockSize),
      static_cast<int>(m_cfg.maxHitsPerTrack));
  ACTS_CUDA_CHECK(cudaDeviceSynchronize());
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_DEBUG("Build tracks...");
  connector.build_tracks();
  ACTS_CUDA_CHECK(cudaDeviceSynchronize());
  ACTS_CUDA_CHECK(cudaGetLastError());

  const int maxHitsPerTrack = connector.cuda_tracks()->max_hits_per_track();
  const int tracksSize = connector.cuda_tracks()->size();
  const int nbTracks = connector.cuda_tracks()->nb_tracks();

  ACTS_DEBUG("maxHitsPerTrack: " << maxHitsPerTrack << ", tracksSize: "
                                 << tracksSize << ", nbTracks: " << nbTracks);

  std::vector<int> nbHits(nbTracks);
  ACTS_CUDA_CHECK(
      cudaMemcpyAsync(nbHits.data(), connector.cuda_tracks()->nb_hits(),
                      nbTracks * sizeof(int), cudaMemcpyDeviceToHost, stream));

  std::vector<int> flatHits(tracksSize);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(
      flatHits.data(), connector.cuda_tracks()->hits(),
      tracksSize * sizeof(int), cudaMemcpyDeviceToHost, stream));

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_CUDA_CHECK(cudaGetLastError());

  std::vector<std::vector<int>> trackCandidates;
  trackCandidates.reserve(nbTracks);

  std::size_t minTrackSize = std::numeric_limits<std::size_t>::max();
  std::size_t maxTrackSize = 0;
  std::size_t avgTrackSize = 0;

  for (int trackIdx = 0; trackIdx < nbTracks; ++trackIdx) {
    const auto* trackBegin = flatHits.data() + trackIdx * maxHitsPerTrack;
    const auto* trackEnd = trackBegin + nbHits[trackIdx];
    trackCandidates.emplace_back(trackBegin, trackEnd);

    // Debug print the first 10 tracks
    if (trackIdx < 10) {
      ACTS_DEBUG("Track " << trackIdx << ": " << [&]() {
        std::ostringstream oss;
        for (const auto& hit : trackCandidates.back()) {
          oss << hit << " ";
        }
        return oss.str();
      }());
    }

    const std::size_t trackSize = trackCandidates.back().size();
    minTrackSize = std::min(minTrackSize, trackSize);
    maxTrackSize = std::max(maxTrackSize, trackSize);
    avgTrackSize += trackSize;
  }

  avgTrackSize /= nbTracks;
  ACTS_DEBUG("Min/Avg/Max track size: " << minTrackSize << "/" << avgTrackSize
                                        << "/" << maxTrackSize);

  return trackCandidates;
}

}  // namespace ActsPlugins
