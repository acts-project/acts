// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>

namespace Acts {

ExaTrkXPipeline::ExaTrkXPipeline(
    std::shared_ptr<GraphConstructionBase> graphConstructor,
    std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
    std::shared_ptr<TrackBuildingBase> trackBuilder,
    std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)),
      m_graphConstructor(graphConstructor),
      m_edgeClassifiers(edgeClassifiers),
      m_trackBuilder(trackBuilder) {
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
    std::vector<int> &spacepointIDs, const ExaTrkXHook &hook,
    ExaTrkXTiming *timing) const {
  try {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto [nodeFeatures, edgeIndex, edgeFeatures] =
        (*m_graphConstructor)(features, spacepointIDs.size(), moduleIds,
                              m_graphConstructor->device());
    auto t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->graphBuildingTime = t1 - t0;
    }

    hook(nodeFeatures, edgeIndex, {});

    std::any edgeScores;
    timing->classifierTimes.clear();

    for (auto edgeClassifier : m_edgeClassifiers) {
      t0 = std::chrono::high_resolution_clock::now();
      auto [newNodeFeatures, newEdgeIndex, newEdgeFeatures, newEdgeScores] =
          (*edgeClassifier)(std::move(nodeFeatures), std::move(edgeIndex),
                            std::move(edgeFeatures), edgeClassifier->device());
      t1 = std::chrono::high_resolution_clock::now();

      if (timing != nullptr) {
        timing->classifierTimes.push_back(t1 - t0);
      }

      nodeFeatures = std::move(newNodeFeatures);
      edgeFeatures = std::move(newEdgeFeatures);
      edgeIndex = std::move(newEdgeIndex);
      edgeScores = std::move(newEdgeScores);

      hook(nodeFeatures, edgeIndex, edgeScores);
    }

    t0 = std::chrono::high_resolution_clock::now();
    auto res = (*m_trackBuilder)(std::move(nodeFeatures), std::move(edgeIndex),
                                 std::move(edgeScores), spacepointIDs,
                                 m_trackBuilder->device());
    t1 = std::chrono::high_resolution_clock::now();

    if (timing != nullptr) {
      timing->trackBuildingTime = t1 - t0;
    }

    return res;
  } catch (Acts::NoEdgesError &) {
    ACTS_WARNING("No egdges left in GNN pipeline, return 0 track candidates");
    if (timing != nullptr) {
      while (timing->classifierTimes.size() < m_edgeClassifiers.size()) {
        timing->classifierTimes.push_back({});
      }
    }
    return {};
  }
}

}  // namespace Acts
