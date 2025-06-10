// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"

#include "Acts/Utilities/Helpers.hpp"

#ifdef ACTS_EXATRKX_WITH_CUDA
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"

namespace {
struct CudaStreamGuard {
  cudaStream_t stream{};
  CudaStreamGuard() { ACTS_CUDA_CHECK(cudaStreamCreate(&stream)); }
  ~CudaStreamGuard() {
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_CUDA_CHECK(cudaStreamDestroy(stream));
  }
};
}  // namespace
#endif

namespace Acts {

ExaTrkXPipeline::ExaTrkXPipeline(
    std::shared_ptr<GraphConstructionBase> graphConstructor,
    std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
    std::shared_ptr<TrackBuildingBase> trackBuilder,
    std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)),
      m_graphConstructor(std::move(graphConstructor)),
      m_edgeClassifiers(std::move(edgeClassifiers)),
      m_trackBuilder(std::move(trackBuilder)) {
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

std::vector<std::vector<int>> ExaTrkXPipeline::run(
    std::vector<float> &features, const std::vector<std::uint64_t> &moduleIds,
    std::vector<int> &spacepointIDs, Acts::Device device,
    const ExaTrkXHook &hook, ExaTrkXTiming *timing) const {
  ExecutionContext ctx;
  ctx.device = device;
#ifdef ACTS_EXATRKX_WITH_CUDA
  std::optional<CudaStreamGuard> streamGuard;
  if (ctx.device.type == Acts::Device::Type::eCUDA) {
    streamGuard.emplace();
    ctx.stream = streamGuard->stream;
  }
#endif

  try {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto tensors =
        (*m_graphConstructor)(features, spacepointIDs.size(), moduleIds, ctx);
    auto t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->graphBuildingTime = t1 - t0;
    }

    hook(tensors, ctx);

    if (timing != nullptr) {
      timing->classifierTimes.clear();
    }

    for (const auto &edgeClassifier : m_edgeClassifiers) {
      t0 = std::chrono::high_resolution_clock::now();
      tensors = (*edgeClassifier)(std::move(tensors), ctx);
      t1 = std::chrono::high_resolution_clock::now();

      if (timing != nullptr) {
        timing->classifierTimes.push_back(t1 - t0);
      }

      hook(tensors, ctx);
    }

    t0 = std::chrono::high_resolution_clock::now();
    auto res = (*m_trackBuilder)(std::move(tensors), spacepointIDs, ctx);
    t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->trackBuildingTime = t1 - t0;
    }

    return res;
  } catch (Acts::NoEdgesError &) {
    ACTS_DEBUG("No edges left in GNN pipeline, return 0 track candidates");
    if (timing != nullptr) {
      while (timing->classifierTimes.size() < m_edgeClassifiers.size()) {
        timing->classifierTimes.push_back({});
      }
    }
    return {};
  }
}

}  // namespace Acts
