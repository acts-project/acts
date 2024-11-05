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

namespace Acts {

class WalkthroughAlgorithm {
 public:
  struct Config {
    float ccScoreCut = 0.01;
    float candidateLowThreshold = 0.1;
    float candidateHighThreshold = 0.6;
    std::size_t minCandidateSize = 3;
  };

  using UndirectedGraph =
      boost::adjacency_list<boost::vecS,         // edge list
                            boost::vecS,         // vertex list
                            boost::undirectedS,  // directedness
                            int,                 // property of vertices
                            float                // property of edges
                            // GraphProperty // graph property
                            >;

  using DirectedGraph = boost::adjacency_list<boost::vecS,       // edge list
                                              boost::vecS,       // vertex list
                                              boost::directedS,  // directedness
                                              int,   // property of vertices
                                              float  // property of edges
                                              // GraphProperty // graph property
                                              >;

  using Vertex = boost::graph_traits<UndirectedGraph>::vertex_descriptor;
  static_assert(std::is_same_v<
                Vertex, boost::graph_traits<DirectedGraph>::vertex_descriptor>);

  using UEdge = boost::graph_traits<UndirectedGraph>::edge_descriptor;
  using DEdge = boost::graph_traits<DirectedGraph>::edge_descriptor;

  WalkthroughAlgorithm(const Config &cfg,
                       std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  template <typename T>
  std::vector<std::vector<int>> operator()(std::span<T> edgeIndexSource,
                                           std::span<T> edgeIndexTarget,
                                           std::span<float> edgeScores,
                                           std::size_t numNodes) const {
    DirectedGraph g(numNodes);

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

    ACTS_VERBOSE("start walkthrough with " << numNodes << " nodes and "
                                           << edgeIndexSource.size()
                                           << " edges");
    connectedComponents(g, finalCandidates, 0);
    ACTS_VERBOSE("Done with walkthrough algorithm, return "
                 << finalCandidates.size() << " candidates");

    return finalCandidates;
  }

 private:
  std::size_t weaklyConnectedComponents(
      const DirectedGraph &g, std::vector<std::size_t> &components) const;
  void connectedComponents(const DirectedGraph &g,
                           std::vector<std::vector<int>> &finalCandidates,
                           int recursion) const;
  void refineSubgraph(const DirectedGraph &graph,
                      std::span<Vertex> subgraphNodes,
                      std::vector<std::vector<int>> &finalCandidates,
                      int recursion) const;
  // void filterEdges(DirectedGraph &graph) const;

  const Acts::Logger &logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts
