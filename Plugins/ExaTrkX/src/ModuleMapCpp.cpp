// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCpp.hpp"

#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"
#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"

#include <graph_creator>
#include <map>
#include <module_map_triplet>

#ifndef ACTS_EXATRKX_CPUONLY
#include <CUDA_graph_creator>
#endif

using namespace torch::indexing;

namespace Acts {

ModuleMapCpp::ModuleMapCpp(const Config &cfg,
                           std::unique_ptr<const Acts::Logger> logger_)
    : m_cfg(cfg), m_logger(std::move(logger_)) {
  if( m_cfg.device == torch::kCPU ) {
    m_graphCreator = std::make_unique<graph_creator<float>>(
        m_cfg.moduleMapPath, 10,
        std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
  } else {
#ifndef ACTS_EXATRKX_CPUONLY
    m_cudaGraphCreator = std::make_unique<CUDA_graph_creator<float>>(
      device.id(), m_cfg.moduleMapPath, 10, std::pair<float, float>{0.f, std::numeric_limits<float>::max()}
    );
#else
    throw std::runtime_error("Cannot use cuda version of GraphModuleMap (CUDA is not enabled in CMake)");
#endif
  }

  if (m_cfg.checkModuleConsistencyPerEvent) {
    for (const auto &[doublet, _] :
         m_graphCreator->get_module_map_triplet().map_doublet()) {
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
}

ModuleMapCpp::~ModuleMapCpp() {}

std::tuple<std::any, std::any, std::any> ModuleMapCpp::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<uint64_t> &moduleIds, torch::Device /*device*/) {
  if (numNodes != moduleIds.size()) {
    throw std::invalid_argument(
        "Module Ids do not match number of graph nodes");
  }

  if (m_cfg.checkModuleConsistencyPerEvent) {
    ACTS_DEBUG("Perform consistency check...");
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

    int hitId = static_cast<int>(i);

    // Needs to be rescaled because ModuleMapGraph expects unscaled features
    float r = hitFeatures[0] * m_cfg.rScale;
    float phi = hitFeatures[1] * m_cfg.phiScale;
    float z = hitFeatures[2] * m_cfg.zScale;

    float x = r * std::cos(phi);
    float y = r * std::sin(phi);

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
  graph<float> graph;
  graph_true<float> graph_true;
  if( m_cfg.device == torch::kCPU ) {
    ACTS_DEBUG("Call CPU graph buildinig");
  std::tie(graph, graph_true) =
      m_graphCreator->build_impl(hitsTree, particlesTree, eventId, print);
  } else {
#ifndef ACTS_EXATRKX_CPUONLY
    CUDA_graph_creator<float>::graph_building_stats stats;
    graph = m_cudaGraphCreator->build_impl(hitsTree, stats);
#else
    throw std::runtime_error("Tried to call CUDA graph buidling which is not enabled");
#endif
  }
  const auto numEdges = boost::num_edges(graph.graph_impl());

  if (numEdges == 0) {
    throw std::runtime_error("no edges");
  }

  ACTS_DEBUG("Got " << numEdges << " edges, put them in a vector...");
  std::vector<int64_t> edgeIndexVector;
  edgeIndexVector.reserve(2 * numEdges);

  constexpr std::size_t numEdgeFeatures = 6;
  std::vector<float> edgeFeatureVector;
  edgeFeatureVector.reserve(numEdgeFeatures * numEdges);

  // TODO I think this is already somewhere in the codebase
  const float pi = static_cast<float>(M_PI);
  auto resetAngle = [pi](float angle) {
    if (angle > pi) {
      return angle - 2.f * pi;
    }
    if (angle < -pi) {
      return angle + 2.f * pi;
    }
    return angle;
  };

  auto [begin, end] = boost::edges(graph.graph_impl());
  for (auto it = begin; it != end; ++it) {
    const auto &edge = *it;

    auto src = graph.graph_impl()[boost::source(edge, graph.graph_impl())];
    auto dst = graph.graph_impl()[boost::target(edge, graph.graph_impl())];

    // Edge index
    assert(src.hit_id() >= 0 && src.hit_id() < static_cast<int>(numNodes));
    assert(dst.hit_id() >= 0 && dst.hit_id() < static_cast<int>(numNodes));

    edgeIndexVector.push_back(src.hit_id());
    edgeIndexVector.push_back(dst.hit_id());

    // Edge features
    // See
    // https://gitlab.cern.ch/gnn4itkteam/acorn/-/blob/dev/acorn/utils/loading_utils.py?ref_type=heads#L288
    const float *srcFeatures = inputValues.data() + src.hit_id() * numFeatures;
    const float *dstFeatures = inputValues.data() + dst.hit_id() * numFeatures;

    const float deltaR = dstFeatures[0] - srcFeatures[0];
    const float deltaPhi =
        resetAngle((dstFeatures[1] - srcFeatures[1]) * m_cfg.phiScale) /
        m_cfg.phiScale;
    const float deltaZ = dstFeatures[2] - srcFeatures[2];
    const float deltaEta = dstFeatures[3] - srcFeatures[3];

    // In acorn, nan gets converted to 0
    float phiSlope = 0.f;
    float rPhiSlope = 0.f;

    if (deltaR != 0) {
      phiSlope = std::clamp(deltaPhi / deltaR, -100.f, 100.f);

      const float avgR = 0.5f * (dstFeatures[0] + srcFeatures[0]);
      rPhiSlope = avgR * phiSlope;
    }

    for (auto f : {deltaR, deltaPhi, deltaZ, deltaEta, phiSlope, rPhiSlope}) {
      edgeFeatureVector.push_back(f);
    }
  }

  // Build final tensors
  ACTS_DEBUG("Construct final tensors...");
  assert(inputValues.size() % numFeatures == 0);
  auto nodeFeatures =
      detail::vectorToTensor2D(inputValues, numFeatures).clone();
  assert(edgeIndexVector.size() % numEdges == 0);
  auto edgeIndex = detail::vectorToTensor2D(edgeIndexVector, 2)
                       .t()
                       .clone()
                       .to(torch::kInt64);

  assert(edgeFeatureVector.size() % numEdgeFeatures == 0);
  auto edgeFeatures =
      detail::vectorToTensor2D(edgeFeatureVector, numEdgeFeatures).clone();

  ACTS_DEBUG("nodeFeatures: " << detail::TensorDetails{nodeFeatures});
  ACTS_DEBUG("edgeIndex: " << detail::TensorDetails{edgeIndex});
  ACTS_DEBUG("edgeFeatures: " << detail::TensorDetails{edgeFeatures});

  ACTS_VERBOSE("Edge vector:\n"
               << (detail::RangePrinter{edgeIndexVector.begin(),
                                        edgeIndexVector.begin() + 10}));
  ACTS_VERBOSE("Edge index slice:\n"
               << edgeIndex.index({Slice(0, 2), Slice(0, 9)}));

  return std::make_tuple(std::move(nodeFeatures), std::move(edgeIndex),
                         std::move(edgeFeatures));
}

}  // namespace Acts
