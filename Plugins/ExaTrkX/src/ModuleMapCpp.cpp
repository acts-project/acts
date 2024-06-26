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
#include <module_map_triplet>

namespace Acts {

ModuleMapCpp::ModuleMapCpp(const Config &cfg,
                           std::unique_ptr<const Acts::Logger> logger_)
    : m_cfg(cfg), m_logger(std::move(logger_)) {
  m_graphCreator = std::make_unique<graph_creator<float>>(
      m_cfg.moduleMapPath, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});

  for (const auto &[doublet, _] :
       m_graphCreator->module_map_triplet().map_doublet()) {
    m_uniqueDoupletModuleIds.push_back(doublet[0]);
    m_uniqueDoupletModuleIds.push_back(doublet[1]);
  }
  std::sort(m_uniqueDoupletModuleIds.begin(), m_uniqueDoupletModuleIds.end());
  auto end = std::unique(m_uniqueDoupletModuleIds.begin(),
                         m_uniqueDoupletModuleIds.end());
  m_uniqueDoupletModuleIds.erase(end, m_uniqueDoupletModuleIds.end());

  ACTS_DEBUG(
      "Unique module IDs in doublets: " << m_uniqueDoupletModuleIds.size());
}

ModuleMapCpp::~ModuleMapCpp() {}

std::tuple<std::any, std::any, std::any> ModuleMapCpp::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<uint64_t> &moduleIds, torch::Device /*device*/) {
  if (numNodes != moduleIds.size()) {
    throw std::invalid_argument(
        "Module Ids do not match number of graph nodes");
  }

  {
    auto uniqueModuleIds = moduleIds;
    std::sort(uniqueModuleIds.begin(), uniqueModuleIds.end());
    auto end = std::unique(uniqueModuleIds.begin(), uniqueModuleIds.end());
    uniqueModuleIds.erase(end, uniqueModuleIds.end());
    ACTS_DEBUG("There are " << uniqueModuleIds.size() << " modules");

    std::vector<uint64_t> moduleIdIntersection;
    moduleIdIntersection.reserve(
        std::max(uniqueModuleIds.size(), m_uniqueDoupletModuleIds.size()));

    std::set_intersection(uniqueModuleIds.begin(), uniqueModuleIds.end(),
                          m_uniqueDoupletModuleIds.begin(),
                          m_uniqueDoupletModuleIds.end(),
                          std::back_inserter(moduleIdIntersection));

    ACTS_DEBUG(
        "Intersection with doublet modules: " << moduleIdIntersection.size());
  }

  const auto numFeatures = inputValues.size() / numNodes;

  hits<float> hitsCollection(false, false);

  ACTS_DEBUG("Start collecting hits...");

  for (auto i = 0ul; i < numNodes; ++i) {
    // TODO Use std::span when we move to C++20
    const float *hitFeatures = inputValues.data() + i * numFeatures;

    uint64_t hitId = i;

    float r = hitFeatures[0];
    float phi = hitFeatures[1];
    float z = hitFeatures[2];

    float x = r*std::cos(phi);
    float y = r*std::sin(phi);

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

  ACTS_DEBUG("Hits tree has " << hitsTree.size()
                              << " hits, now build graph...");
  bool print = logger().level() == Acts::Logging::VERBOSE;
  auto [graph, _] =
      m_graphCreator->build_impl(hitsTree, particlesTree, eventId, print);
  const auto numEdges = boost::num_edges(graph.graph_impl());

  if (numEdges == 0) {
    throw std::runtime_error("no edges");
  }

  ACTS_DEBUG("Got " << numEdges << " edges, put them in a vector...");
  std::vector<int32_t> edgeIndexVector;
  edgeIndexVector.reserve(2 * numEdges);

  std::vector<float> edgeFeatureVector;

  auto [begin, end] = boost::edges(graph.graph_impl());
  for (auto it = begin; it != end; ++it) {
    const auto &edge = *it;

    auto src = graph.graph_impl()[boost::source(edge, graph.graph_impl())];
    auto tgt = graph.graph_impl()[boost::target(edge, graph.graph_impl())];

    // Edge index
    assert(src.hit_id() >= 0 && src.hit_id() < static_cast<int>(numNodes));
    assert(tgt.hit_id() >= 0 && tgt.hit_id() < static_cast<int>(numNodes));

    edgeIndexVector.push_back(src.hit_id());
    edgeIndexVector.push_back(tgt.hit_id());

    // Edge features
    edgeFeatureVector.push_back(graph.graph_impl()[edge].dr());
    edgeFeatureVector.push_back(graph.graph_impl()[edge].dPhi());
    edgeFeatureVector.push_back(graph.graph_impl()[edge].dz());
    edgeFeatureVector.push_back(graph.graph_impl()[edge].dEta());

    const auto deltaR = tgt.r() - src.r();
    const auto deltaPhi = tgt.phi() - src.phi();
    const auto phiSlope = deltaPhi / deltaR;
    edgeFeatureVector.push_back(phiSlope);

    const auto deltaRPhi = tgt.r()*tgt.phi() - src.r()*src.phi();
    const float rPhiSlope = deltaRPhi / deltaR;
    edgeFeatureVector.push_back(rPhiSlope);
  }

  // Build final tensors
  ACTS_DEBUG("Construct final tensors...");
  auto nodeFeatures = detail::vectorToTensor2D(inputValues, numFeatures);
  auto edgeIndex = detail::vectorToTensor2D(edgeIndexVector, numEdges);

  constexpr std::size_t numEdgeFeatures = 6;
  auto edgeFeatures =
      detail::vectorToTensor2D(edgeFeatureVector, numEdgeFeatures);

  ACTS_DEBUG("nodeFeatures: " << nodeFeatures.sizes()
                              << ", edgeIndex: " << edgeIndex.sizes()
                              << ", edgeFeatures: " << edgeFeatures.sizes());

  return std::make_tuple(std::move(nodeFeatures), std::move(edgeIndex),
                         std::move(edgeFeatures));
}

}  // namespace Acts
