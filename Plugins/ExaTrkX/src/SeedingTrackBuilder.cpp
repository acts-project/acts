// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/SeedingTrackBuilder.hpp"

#include "Acts/Utilities/Zip.hpp"

#include <map>
#include <ranges>
#include <span>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <torch/torch.h>

namespace stdr = std::ranges;
using namespace torch::indexing;

template <typename directedness_t, typename weight_t>
using Graph = boost::adjacency_list<boost::vecS,         // edge list
                                    boost::vecS,         // vertex list
                                    directedness_t,      // directedness
                                    boost::no_property,  // property of vertices
                                    weight_t             // property of edges
                                    >;

namespace Acts {

std::vector<std::vector<int>> SeedingTrackBuilder::operator()(
    std::any nodes, std::any edges, std::any weights,
    std::vector<int> &spacepointIDs, torch::Device) {
  ACTS_DEBUG("Start track building");
  const auto nodeTensor = std::any_cast<torch::Tensor>(nodes).to(torch::kCPU);
  const auto rTensor =
      nodeTensor.index({None, m_cfg.rColumn}).contiguous().to(torch::kFloat);
  const auto zTensor =
      nodeTensor.index({None, m_cfg.zColumn}).contiguous().to(torch::kFloat);

  const auto theta = torch::atan(rTensor / zTensor);
  const auto eta = -torch::log(torch::tan(0.5 * theta));

  const auto edgeTensor = std::any_cast<torch::Tensor>(edges).to(torch::kCPU);
  const auto edgeWeightTensor =
      std::any_cast<torch::Tensor>(weights).to(torch::kCPU);

  assert(edgeTensor.size(0) == 2);
  assert(edgeTensor.size(1) == edgeWeightTensor.size(0));

  const auto numSpacepoints = spacepointIDs.size();
  const auto numEdges = static_cast<std::size_t>(edgeWeightTensor.size(0));

  if (numEdges == 0) {
    ACTS_WARNING("No edges remained after edge classification");
    return {};
  }

  using vertex_t = int;
  using weight_t = float;

  std::span<float> rValues(rTensor.data_ptr<float>(), numSpacepoints);
  std::span<float> etaValues(eta.data_ptr<float>(), numSpacepoints);

  std::span<vertex_t> rowIndices(edgeTensor.data_ptr<vertex_t>(), numEdges);
  std::span<vertex_t> colIndices(edgeTensor.data_ptr<vertex_t>() + numEdges,
                                 numEdges);
  std::span<weight_t> edgeWeights(edgeWeightTensor.data_ptr<float>(), numEdges);

  using UndirectedGraph = Graph<boost::undirectedS, weight_t>;
  using DirectedGraph = Graph<boost::bidirectionalS, boost::no_property>;

  UndirectedGraph graph(numSpacepoints);

  for (const auto [row, col, weight] :
       Acts::zip(rowIndices, colIndices, edgeWeights)) {
    auto avgEta = 0.5 * (etaValues[row] + etaValues[col]);
    float scoreCut = m_cfg.baseScoreCut;

    for (const auto &[etaMin, etaMax, regionScoreCut] : m_cfg.scoreCutRegions) {
      if (avgEta >= etaMin && avgEta < etaMax) {
        scoreCut = regionScoreCut;
      }
    }

    if (weight > scoreCut) {
      boost::add_edge(row, col, weight, graph);
    }
  }

  std::vector<vertex_t> labels(numSpacepoints);
  auto numLabels = boost::connected_components(graph, &labels[0]);

  std::vector<DirectedGraph> componentGraphs(numLabels);
  const auto [eBegin, eEnd] = boost::edges(graph);

  // Because the components are connected, both vertices for each edge must have
  // the same label by definition
  for (auto edgeIt = eBegin; edgeIt != eEnd; ++edgeIt) {
    vertex_t src = boost::source(*edgeIt, graph);
    vertex_t tgt = boost::target(*edgeIt, graph);
    assert(labels.at(src) == labels.at(tgt));

    // We make the graph directed dependent on radius r (should be more reliable
    // then sqrt(r^2 + z^2))
    if (rValues[src] > rValues[tgt]) {
      std::swap(src, tgt);
    }

    const auto label = labels.at(src);
    auto &componentGraph = componentGraphs.at(label);

    boost::add_edge(src, tgt, componentGraph);
  }

  std::vector<std::vector<vertex_t>> trackCandidates;

  // TODO copying the track candidates each recursion is expensive, but
  // optimization should come at a later stage when it actually works...
  auto buildTrackRecursive = [&](const auto &self,
                                 std::vector<vertex_t> trackCandidate,
                                 const auto &subgraph) -> void {
    const auto src = trackCandidate.back();

    std::size_t nextVtcs = 0;
    auto [begin, end] = boost::out_edges(src, subgraph);
    for (auto edge = begin; edge != end; ++edge) {
      const auto tgt = boost::target(*edge, subgraph);
      auto newCand = trackCandidate;
      newCand.push_back(tgt);

      self(self, newCand, subgraph);
      nextVtcs++;
    }

    if (nextVtcs == 0) {
      if (trackCandidate.size() > m_cfg.minSeedSize) {
        ACTS_VERBOSE("Built track candidate with " << trackCandidate.size()
                                                   << " hits");
        trackCandidates.push_back(trackCandidate);
      } else {
        ACTS_VERBOSE("Discard track candidate with " << trackCandidate.size()
                                                     << " hits");
      }
    }
  };

  std::size_t nMoreThanOneCandidatePerLabel = 0;
  for (auto label = 0ul; label < componentGraphs.size(); ++label) {
    const auto &subgraph = componentGraphs[label];

    const auto nCandidatesBefore = trackCandidates.size();
    std::size_t nStartVertices = 0;

    auto [vbegin, vend] = boost::vertices(subgraph);
    for (auto vtxIt = vbegin; vtxIt != vend; ++vtxIt) {
      vertex_t vtx = *vtxIt;
      if (boost::in_degree(vtx, subgraph) == 0) {
        nStartVertices++;
        buildTrackRecursive(buildTrackRecursive, {vtx}, subgraph);
      }
    }

    const auto nCandidatesFound = trackCandidates.size() - nCandidatesBefore;
    if (nCandidatesFound > 0) {
      nMoreThanOneCandidatePerLabel++;
    }
    ACTS_VERBOSE("Label " << label << ": For " << nStartVertices << " found "
                          << trackCandidates.size() - nCandidatesBefore
                          << " track candidates");
  }

  ACTS_DEBUG("Connected-component-labels with more than one candidate: "
             << nMoreThanOneCandidatePerLabel << " / " << numLabels);
  ACTS_DEBUG("Found track candidates: " << trackCandidates.size());
}

}  // namespace Acts
