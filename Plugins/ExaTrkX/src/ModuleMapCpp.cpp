// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCpp.hpp"

#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"

#include <graph_creator>

namespace Acts {

ModuleMapCpp::ModuleMapCpp(const Config &cfg,
                           std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  m_graphCreator = std::make_unique<graph_creator<float>>(
      m_cfg.moduleMapPath, 10, std::pair<float, float>{0.f, 10.f}, std::nullopt,
      std::nullopt, false, false);
}

ModuleMapCpp::~ModuleMapCpp() {}

std::tuple<std::any, std::any> ModuleMapCpp::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<uint64_t> &moduleIds, int /*deviceHint*/) {
  if (numNodes != moduleIds.size()) {
    throw std::invalid_argument(
        "Module Ids do not match number of graph nodes");
  }

  const auto numFeatures = inputValues.size() / numNodes;

  hits<float> hitsCollection(false, false);

  for (auto i = 0ul; i < numNodes; ++i) {
    // TODO Use std::span when we move to C++20
    const float *hitFeatures = inputValues.data() + i * numFeatures;

    uint64_t hitId = i;

    // TODO this is currently wrong!!! we do not have x, y, z here usually
    // Somehow we need to provide a mecanism to give x, y, z to the algorithm
    float x = hitFeatures[0];
    float y = hitFeatures[1];
    float z = hitFeatures[2];

    uint64_t particleId = 0;  // We do not know
    uint64_t moduleId = moduleIds[i];
    std::string hardware = "";  // now hardware
    int barrelEndcap = 0;       // unclear, is this a flag???
    uint64_t particleID1 = 0;   // unclear
    uint64_t particleID2 = 0;   // unclear

    hit<float> hit(hitId, x, y, z, particleId, moduleId, hardware, barrelEndcap,
                   particleID1, particleID2);

    hitsCollection += hit;
  }

  TTree_hits<float> hitsTree = hitsCollection;
  TTree_particles<float> particlesTree;  // dummy, not needed currently
  std::string eventId = "no-id";

  auto [graph, _] =
      m_graphCreator->build_impl(hitsTree, particlesTree, eventId, false);
  const auto numEdges = boost::num_edges(graph.graph_impl());

  std::vector<int32_t> edgeVector;
  edgeVector.reserve(2 * numEdges);

  auto [begin, end] = boost::edges(graph.graph_impl());
  for (auto it = begin; it != end; ++it) {
    const auto &edge = *it;

    auto source = boost::source(edge, graph.graph_impl());
    auto target = boost::target(edge, graph.graph_impl());

    edgeVector.push_back(graph.graph_impl()[source].hit_id());
    edgeVector.push_back(graph.graph_impl()[target].hit_id());
  }

  // Build final tensors
  auto featureTensor = detail::vectorToTensor2D(inputValues, numFeatures);
  auto edgeTensor = detail::vectorToTensor2D(edgeVector, numEdges);

  return std::make_tuple(std::move(featureTensor), std::move(edgeTensor));
}

}  // namespace Acts
