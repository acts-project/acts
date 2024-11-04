// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Zip.hpp>

#include <span>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/topological_sort.hpp>

#include "Stages.hpp"

namespace Acts {

class WalkthroughAlgorithm {
 public:
  struct Config {
    float ccScoreCut = 0.01;
    float highScoreCut = 0.6;
    std::size_t minCandidateSize = 3;
  };

  using Graph = boost::adjacency_list<boost::vecS,       // edge list
                                      boost::vecS,       // vertex list
                                      boost::directedS,  // directedness
                                      int,               // property of vertices
                                      float              // property of edges
                                      >;

  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = boost::graph_traits<Graph>::edge_descriptor;

  WalkthroughAlgorithm(const Config &cfg,
                       std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  template <typename T>
  std::vector<std::vector<int>> operator()(std::span<T> edgeIndexSource,
                                           std::span<T> edgeIndexTarget,
                                           std::span<float> edgeScores,
                                           std::size_t numNodes) const {
    Graph g(numNodes);

    for (const auto [source, target, weight] :
         Acts::zip(edgeIndexSource, edgeIndexTarget, edgeScores)) {
      if (weight > m_cfg.ccScoreCut) {
        boost::add_edge(source, target, weight, g);
      }
    }

    for (Vertex v = 0; v < numNodes; ++v) {
      g[v] = v;
    }

    std::vector<std::vector<int>> finalCandidates;

    connectedComponents(g, finalCandidates);

    return finalCandidates;
  }

 private:
  void connectedComponents(
      const Graph &g, std::vector<std::vector<int>> &finalCandidates) const;
  void refineSubgraph(const Graph &graph, std::span<Vertex> subgraphNodes,
                      std::vector<std::vector<int>> &finalCandidates) const;
  void filterEdges(Graph &graph) const;

  const Acts::Logger &logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

class BoostWalkthrough final : public Acts::TrackBuildingBase {
 public:
  BoostWalkthrough(std::unique_ptr<const Logger> logger)
      : m_logger(std::move(logger)), m_device(torch::Device(torch::kCPU)) {}

  std::vector<std::vector<int>> operator()(
      std::any nodes, std::any edges, std::any edge_weights,
      std::vector<int> &spacepointIDs,
      torch::Device device = torch::Device(torch::kCPU)) override;
  torch::Device device() const override { return m_device; };

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  torch::Device m_device;
  const auto &logger() const { return *m_logger; }
};

}  // namespace Acts
