// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"

#include "Acts/Utilities/Zip.hpp"

#include <algorithm>
#include <map>
#include <span>

#include <boost/beast/core/span.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <torch/torch.h>

namespace bc = boost::container;

namespace {
template <typename vertex_t, typename weight_t>
auto weaklyConnectedComponents(vertex_t numNodes,
                               std::span<vertex_t>& srcIndices,
                               std::span<vertex_t>& dstIndices,
                               std::span<weight_t>& edgeWeights,
                               std::vector<vertex_t>& trackLabels) {
  using Graph =
      boost::adjacency_list<boost::vecS,         // edge list
                            boost::vecS,         // vertex list
                            boost::undirectedS,  // directedness
                            boost::no_property,  // property of vertices
                            weight_t             // property of edges
                            >;

  Graph g(numNodes);

  for (const auto [row, col, weight] :
       Acts::zip(srcIndices, dstIndices, edgeWeights)) {
    boost::add_edge(row, col, weight, g);
  }

  return boost::connected_components(g, &trackLabels[0]);
}
}  // namespace

namespace Acts {

std::vector<std::vector<int>> BoostTrackBuilding::operator()(
    std::any /*nodes*/, std::any edges, std::any weights,
    std::vector<int>& spacepointIDs, const ExecutionContext& execContext) {
  ACTS_DEBUG("Start track building");
  const auto edgeTensor = std::any_cast<torch::Tensor>(edges).to(torch::kCPU);
  const auto edgeWeightTensor =
      std::any_cast<torch::Tensor>(weights).to(torch::kCPU);

  assert(edgeTensor.size(0) == 2);
  assert(edgeTensor.size(1) == edgeWeightTensor.size(0));

  const auto numSpacepoints = spacepointIDs.size();
  const auto numEdgesInitial =
      static_cast<std::size_t>(edgeWeightTensor.size(0));

  if (numEdgesInitial == 0) {
    ACTS_WARNING("No edges remained after edge classification");
    return {};
  }

  using vertex_t = std::int64_t;
  using weight_t = float;

  std::span<vertex_t> srcIndices(edgeTensor.data_ptr<vertex_t>(),
                                 numEdgesInitial);
  std::span<vertex_t> dstIndices(
      edgeTensor.data_ptr<vertex_t>() + numEdgesInitial, numEdgesInitial);
  std::span<weight_t> edgeWeights(edgeWeightTensor.data_ptr<float>(),
                                  numEdgesInitial);

  if (m_cfg.doWalkthrough) {
    ACTS_DEBUG("Do walkthrough algorithm");
    return m_walkthrough(srcIndices, dstIndices, edgeWeights, numSpacepoints);
  }

  // In case we do junction removal, these vectors hold the resulting new edges
  std::vector<vertex_t> srcIndicesJR, dstIndicesJR;
  std::vector<weight_t> edgeWeightsJR;

  if (m_cfg.doJunctionRemoval) {
    ACTS_DEBUG("Do junction removal before connected components");

    // Collect the indices into the edge lists for incoming and ougoing edges
    // per spacepoint
    std::vector<bc::small_vector<std::size_t, 4>> incomingEdges(numSpacepoints);
    std::vector<bc::small_vector<std::size_t, 4>> outgoingEdges(numSpacepoints);

    std::size_t maxOutEdges{}, maxInEdges{};
    for (std::size_t i = 0; i < numEdgesInitial; ++i) {
      auto src = srcIndices[i];
      auto dst = dstIndices[i];

      outgoingEdges.at(src).push_back(i);
      incomingEdges.at(dst).push_back(i);

      maxOutEdges = std::max(maxOutEdges, outgoingEdges.at(src).size());
      maxInEdges = std::max(maxInEdges, incomingEdges.at(dst).size());
    }

    ACTS_DEBUG("JR max incoming edges: "
               << maxInEdges << ", max outgoing edges: " << maxOutEdges);

    // In case more than 1 incoming/outgoing edges, only keep the best
    // Use a boolean mask to store the result
    std::vector<bool> edgeMask(numEdgesInitial, true);
    for (std::size_t i = 0; i < numSpacepoints; ++i) {
      for (const auto& ioEdges : {incomingEdges.at(i), outgoingEdges.at(i)}) {
        if (ioEdges.size() < 2) {
          continue;
        }
        auto maxEdge = *std::ranges::max_element(
            ioEdges, {}, [&](auto i) { return edgeWeights[i]; });
        for (auto e : ioEdges) {
          if (e != maxEdge) {
            edgeMask.at(e) = false;
          }
        }
      }
    }

    // Copy the remaining edges into new vectors
    srcIndicesJR.reserve(numEdgesInitial);
    dstIndicesJR.reserve(numEdgesInitial);

    for (auto [src, dst, weight, keep] :
         Acts::zip(srcIndices, dstIndices, edgeWeights, edgeMask)) {
      if (keep) {
        srcIndicesJR.push_back(src);
        dstIndicesJR.push_back(dst);
        edgeWeightsJR.push_back(weight);
      }
    }

    ACTS_DEBUG("JR removed " << numEdgesInitial - srcIndicesJR.size() << " / "
                             << numEdgesInitial << " edges");

    // Replace the spans so that they point to the new edge list
    srcIndices = std::span<vertex_t>(srcIndicesJR.data(), srcIndicesJR.size());
    dstIndices = std::span<vertex_t>(dstIndicesJR.data(), dstIndicesJR.size());
    edgeWeights =
        std::span<weight_t>(edgeWeightsJR.data(), edgeWeightsJR.size());
  }

  ACTS_DEBUG("Do simple connected components");
  std::vector<vertex_t> trackLabels(numSpacepoints);

  auto numberLabels = weaklyConnectedComponents<vertex_t, weight_t>(
      numSpacepoints, srcIndices, dstIndices, edgeWeights, trackLabels);

  ACTS_VERBOSE("Number of track labels: " << trackLabels.size());
  ACTS_VERBOSE("Number of unique track labels: " << [&]() {
    std::vector<vertex_t> sorted(trackLabels);
    std::ranges::sort(sorted);
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    return sorted.size();
  }());

  if (trackLabels.size() == 0) {
    return {};
  }

  std::vector<std::vector<int>> trackCandidates(numberLabels);

  for (const auto [label, id] : Acts::zip(trackLabels, spacepointIDs)) {
    trackCandidates[label].push_back(id);
  }

  return trackCandidates;
}

}  // namespace Acts
