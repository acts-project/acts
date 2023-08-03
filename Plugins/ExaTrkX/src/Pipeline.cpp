// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/Pipeline.hpp"

namespace Acts {

Pipeline::Pipeline(
    std::shared_ptr<GraphConstructionBase> graphConstructor,
    std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
    std::shared_ptr<TrackBuildingBase> trackBuilder, const Acts::Logger &logger)
    : m_logger(logger.cloneWithSuffix("Pipeline")),
      m_graphConstructor(graphConstructor),
      m_edgeClassifiers(edgeClassifiers),
      m_trackBuilder(trackBuilder) {
  if (!m_graphConstructor) {
    throw std::invalid_argument("Missing graph construction module");
  }
  if (!m_trackBuilder) {
    throw std::invalid_argument("Missing track building module");
  }
  if (m_edgeClassifiers.empty() or
      not std::all_of(m_edgeClassifiers.begin(), m_edgeClassifiers.end(),
                      [](const auto &a) { return static_cast<bool>(a); })) {
    throw std::invalid_argument("Missing graph construction module");
  }
}

std::vector<std::vector<int>> Pipeline::run(
    boost::multi_array<float, 2> &features, std::vector<int> &spacepointIDs,
    const PipelineHook &hook) const {
  auto [nodes, edges] = (*m_graphConstructor)(features);

  hook(nodes, edges, logger());

  std::any edge_weights;

  for (auto edgeClassifier : m_edgeClassifiers) {
    auto [newNodes, newEdges, newWeights] =
        (*edgeClassifier)(std::move(nodes), std::move(edges));
    nodes = std::move(newNodes);
    edges = std::move(newEdges);
    edge_weights = std::move(newWeights);

    hook(nodes, edges, logger());
  }

  return (*m_trackBuilder)(std::move(nodes), std::move(edges),
                           std::move(edge_weights), spacepointIDs);
}

}  // namespace Acts
