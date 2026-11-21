// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/GnnPipeline.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsPlugins/Gnn/detail/NvtxUtils.hpp"

#include <optional>

#ifdef ACTS_GNN_WITH_CUDA
#include "ActsPlugins/Gnn/detail/async_device_memory_resource.hpp"

#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>
#endif

using namespace Acts;

namespace ActsPlugins {

GnnPipeline::GnnPipeline(
    std::shared_ptr<GraphConstructionBase> graphConstructor,
    std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
    std::shared_ptr<TrackBuildingBase> trackBuilder,
    std::unique_ptr<const Logger> logger, bool shrinkNodes)
    : m_logger(std::move(logger)),
      m_graphConstructor(std::move(graphConstructor)),
      m_edgeClassifiers(std::move(edgeClassifiers)),
      m_trackBuilder(std::move(trackBuilder)),
      m_shrinkNodes(shrinkNodes) {
  if (!m_graphConstructor) {
    throw std::invalid_argument("Missing graph construction module");
  }
  if (!m_trackBuilder) {
    throw std::invalid_argument("Missing track building module");
  }
  if (m_edgeClassifiers.empty() ||
      rangeContainsValue(m_edgeClassifiers, nullptr)) {
    throw std::invalid_argument("Missing graph construction module");
  }
}

std::vector<std::vector<int>> GnnPipeline::run(
    std::vector<float> &features, const std::vector<std::uint64_t> &moduleIds,
    std::vector<int> &spacePointIds, Device device, const GnnHook &hook,
    GnnTiming *timing) const {
  ExecutionContext ctx;
  ctx.device = device;
#ifdef ACTS_GNN_WITH_CUDA
  // A new stream, memory resource and copy object for every run, so that runs
  // on several threads stay independent
  std::optional<vecmem::cuda::stream_wrapper> stream;
  std::optional<detail::async_device_memory_resource> deviceMr;
  std::optional<vecmem::cuda::async_copy> copy;
  if (ctx.device.type == Device::Type::eCUDA) {
    stream.emplace();
    deviceMr.emplace(static_cast<cudaStream_t>(stream->stream()));
    copy.emplace(*stream);
    ctx.stream = static_cast<cudaStream_t>(stream->stream());
    ctx.memoryResource = &*deviceMr;
    ctx.copy = &*copy;
  }
#endif

  try {
    auto t0 = std::chrono::high_resolution_clock::now();
    ACTS_NVTX_START(graph_construction);
    auto tensors =
        (*m_graphConstructor)(features, spacePointIds.size(), moduleIds, ctx);
    ACTS_NVTX_STOP(graph_construction);
    auto t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->graphBuildingTime = t1 - t0;
    }

    hook(tensors, ctx);

    if (timing != nullptr) {
      timing->classifierTimes.clear();
    }

    for (const auto &edgeClassifier : m_edgeClassifiers) {
      if (m_shrinkNodes) {
        const auto nodesBefore = tensors.nodeFeatures.shape()[0];
        tensors = removeUnusedNodes(std::move(tensors), spacePointIds, ctx);
        ACTS_DEBUG("Removed unused nodes: " << nodesBefore << " -> "
                                            << tensors.nodeFeatures.shape()[0]);
      }

      t0 = std::chrono::high_resolution_clock::now();
      ACTS_NVTX_START(edge_classifier);
      tensors = (*edgeClassifier)(std::move(tensors), ctx);
      ACTS_NVTX_STOP(edge_classifier);
      t1 = std::chrono::high_resolution_clock::now();

      if (timing != nullptr) {
        timing->classifierTimes.push_back(t1 - t0);
      }

      hook(tensors, ctx);
    }

    t0 = std::chrono::high_resolution_clock::now();
    ACTS_NVTX_START(track_building);
    auto res = (*m_trackBuilder)(std::move(tensors), spacePointIds, ctx);
    ACTS_NVTX_STOP(track_building);
    t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->trackBuildingTime = t1 - t0;
    }

    return res;
  } catch (NoEdgesError &) {
    ACTS_DEBUG("No edges left in GNN pipeline, return 0 track candidates");
    if (timing != nullptr) {
      while (timing->classifierTimes.size() < m_edgeClassifiers.size()) {
        timing->classifierTimes.push_back({});
      }
    }
    return {};
  }
}

}  // namespace ActsPlugins
