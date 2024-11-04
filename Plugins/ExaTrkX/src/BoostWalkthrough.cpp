// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/BoostWalkthrough.hpp"

#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Zip.hpp>

#include <span>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/topological_sort.hpp>

namespace {

template <typename Graph>
struct EdgePredicate {
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

  const Graph *graph{};
  const std::span<Vertex> *vertices{};

  template <typename Edge>
  bool operator()(const Edge &e) const {
    return std::ranges::find(*vertices, boost::source(e, *graph)) !=
           vertices->end();
  }
};

template <typename Graph>
struct VertexPredicate {
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

  const std::span<Vertex> *vertices{};

  template <typename Vertex>
  bool operator()(const Vertex &v) const {
    return std::ranges::find(*vertices, v) != vertices->end();
  }
};

template <typename Graph>
using FilteredGraph =
    boost::filtered_graph<Graph, EdgePredicate<Graph>, VertexPredicate<Graph>>;

template <typename Graph,
          typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
Graph extractSubgraph(const Graph &graph, std::span<Vertex> subgraphVertices) {
  EdgePredicate<Graph> edgePred{&graph, &subgraphVertices};
  VertexPredicate<Graph> vertPred{&subgraphVertices};

  FilteredGraph<Graph> filteredGraph(graph, edgePred, vertPred);
  Graph subgraph;
  boost::copy_graph(filteredGraph, subgraph);

  return subgraph;
}

template <typename Graph,
          typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::vector<int> transformCandidate(const std::vector<Vertex> &vec,
                                    const Graph &graph) {
  std::vector<int> ret(vec.size());
  std::ranges::transform(vec, ret.begin(), [&](Vertex v) { return graph[v]; });
  return ret;
}

}  // namespace

namespace Acts {

void WalkthroughAlgorithm::connectedComponents(
    const Graph &g, std::vector<std::vector<int>> &finalCandidates) const {
  // 2) Connected components
  std::vector<Vertex> components(boost::num_vertices(g));
  auto numComponents = boost::connected_components(g, components.data());
  ACTS_VERBOSE("Connected component resulted in " << numComponents
                                                  << " components");

  // 3) Search for chain like components. Do this by going once through all
  // vertices and store for each chain the number of chain-like nodes
  // There is a bit reduncancy in this struct, but build it for debugging
  struct ChainProperties {
    int nOneOutgoing = 0;
    int nZeroOutgoing = 0;
    int nOther = 0;
  };

  std::vector<ChainProperties> chainProperties(numComponents);
  std::vector<std::vector<Vertex>> componentNodes(numComponents);

  for (Vertex v = 0; v < boost::num_vertices(g); ++v) {
    auto component = components.at(v);
    auto outDegree = boost::out_degree(v, g);
    if (outDegree == 0) {
      chainProperties.at(component).nZeroOutgoing++;
    } else if (outDegree == 1) {
      chainProperties.at(component).nOneOutgoing++;
    } else {
      chainProperties.at(component).nOther++;
    }

    componentNodes.at(component).push_back(v);
  }

  for (auto c = 0ul; c < numComponents; ++c) {
    const auto &props = chainProperties.at(c);
    if (componentNodes.at(c).size() < m_cfg.minCandidateSize) {
      continue;
    }
    if (props.nOther == 0 && props.nZeroOutgoing == 1) {
      ACTS_VERBOSE("Found chain like component with length "
                   << componentNodes.at(c).size());
      finalCandidates.emplace_back(transformCandidate(componentNodes.at(c), g));
    } else {
      refineSubgraph(g, componentNodes.at(c), finalCandidates);
    }
  }
}

void WalkthroughAlgorithm::refineSubgraph(
    const Graph &graph, std::span<Vertex> subgraphNodes,
    std::vector<std::vector<int>> &finalCandidates) const {
  auto subgraph = extractSubgraph(graph, subgraphNodes);

  // Make a graph that only retains the edge with the largest score, or edges
  // with a score above highScoreCut
  filterEdges(subgraph);

  ACTS_VERBOSE("Refine subgraph with "
               << boost::num_vertices(subgraph) << " vertices and "
               << boost::num_edges(subgraph) << " edges");

  // Do a topological sort to find a source node
  std::vector<Vertex> topologicalSortedReversed(subgraphNodes.size());
  try {
    boost::topological_sort(subgraph, topologicalSortedReversed.begin());
  } catch (boost::not_a_dag &) {
    ACTS_ERROR("Encountered acyclic graph, no candidates are produced!");
    return;
  }

  Vertex sourceNode = topologicalSortedReversed.back();
  ACTS_VERBOSE("Source node of path building is " << sourceNode);

  std::vector<std::pair<bool, std::vector<Vertex>>> paths;
  paths.push_back({false, {sourceNode}});

  std::size_t finishedPaths = 0;
  while (finishedPaths < paths.size()) {
    std::vector<std::pair<bool, std::vector<Vertex>>> newPaths;

    for (auto &[finished, path] : paths) {
      if (finished) {
        continue;
      }

      auto [begin, end] = boost::out_edges(path.back(), subgraph);
      auto outDegree = std::distance(begin, end);

      if (outDegree == 0) {
        finishedPaths++;
        finished = true;
        continue;
      }

      // Extend the current path with the first outgoing edge
      path.push_back(boost::target(*begin, subgraph));

      // Loop over the remaining edges and create new path candidates
      for (auto it = std::next(begin); it != end; ++it) {
        std::vector<Vertex> newPath(path.begin(), std::prev(path.end()));
        newPath.push_back(boost::target(*it, subgraph));

        newPaths.emplace_back(false, std::move(newPath));
      }
    }

    // Add new paths to path collection
    for (auto &newPath : newPaths) {
      paths.emplace_back(std::move(newPath));
    }
    ACTS_VERBOSE("- Found " << newPaths.size() << " new paths, now have "
                            << paths.size() << " paths");
  }

  for (const auto &path : paths) {
    ACTS_VERBOSE("path: " << [&]() {
      std::stringstream ss;
      for (const auto &v : path.second) {
        ss << subgraph[v] << " ";
      }
      return ss.str();
    }());
  }

  // Add the longest candidate to the final list
  auto longestPathIt = std::ranges::max_element(
      paths, std::less{}, [](const auto &p) { return p.second.size(); });
  auto &finalPath = longestPathIt->second;
  ACTS_VERBOSE("Found " << paths.size() << " paths, longest path has length "
                        << finalPath.size());

  finalCandidates.emplace_back(transformCandidate(finalPath, subgraph));

  if (finalPath.size() == boost::num_vertices(subgraph)) {
    ACTS_VERBOSE("Done with subgraph, no nodes left");
    return;
  }

  std::vector<Vertex> remainingNodes(boost::num_vertices(subgraph) -
                                     finalPath.size());
  auto [vbegin, vend] = boost::vertices(subgraph);
  std::set_difference(vbegin, vend, finalPath.begin(), finalPath.end(),
                      remainingNodes.begin());

  // Do connected components and recurse
  auto nextGraph =
      extractSubgraph(subgraph, {remainingNodes.begin(), remainingNodes.end()});
  connectedComponents(nextGraph, finalCandidates);
}

void Acts::WalkthroughAlgorithm::filterEdges(Graph &graph) const {
  auto [vbegin, vend] = boost::vertices(graph);
  for (auto vit = vbegin; vit != vend; ++vit) {
    auto [begin, end] = boost::out_edges(*vit, graph);
    float maxScore = 0.f;

    for (auto it = begin; it != end; ++it) {
      maxScore = std::max(graph[*it], maxScore);
    }

    for (auto it = begin; it != end; ++it) {
      if (!(graph[*it] == maxScore || graph[*it] > m_cfg.highScoreCut)) {
        graph[*it] = 0.f;
      }
    }
  }

  boost::remove_edge_if([&](auto e) { return graph[e] == 0.f; }, graph);
}

}  // namespace Acts
