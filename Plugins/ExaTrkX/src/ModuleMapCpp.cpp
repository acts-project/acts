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

std::tuple<std::any, std::any> ModuleMapCpp::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<uint64_t> &moduleIds, int /*deviceHint*/) {
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

    // Make sure that the features are provided as x,y,z and not r,phi,z
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
  std::vector<int32_t> edgeVector;
  edgeVector.reserve(2 * numEdges);

  auto [begin, end] = boost::edges(graph.graph_impl());
  for (auto it = begin; it != end; ++it) {
    const auto &edge = *it;

    auto source = boost::source(edge, graph.graph_impl());
    auto target = boost::target(edge, graph.graph_impl());

    auto sid = graph.graph_impl()[source].hit_id();
    auto tid = graph.graph_impl()[target].hit_id();

    assert(sid >= 0 && sid < static_cast<int>(numNodes));
    assert(tid >= 0 && tid < static_cast<int>(numNodes));

    edgeVector.push_back(sid);
    edgeVector.push_back(tid);
  }

  // Build final tensors
  ACTS_DEBUG("Construct final tensors...");
  auto featureTensor = detail::vectorToTensor2D(inputValues, numFeatures);
  auto edgeTensor = detail::vectorToTensor2D(edgeVector, numEdges);

  ACTS_DEBUG(
      "All edges < numNodes: "
      << std::boolalpha
      << (edgeTensor < static_cast<int32_t>(numNodes)).all().item<bool>());

  ACTS_DEBUG("featureTensor: " << featureTensor.sizes()
                               << ", edgeTensor: " << edgeTensor.sizes());
  return std::make_tuple(std::move(featureTensor), std::move(edgeTensor));
}

}  // namespace Acts
