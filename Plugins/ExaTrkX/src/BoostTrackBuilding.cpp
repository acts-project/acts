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

#include <boost/beast/core/span.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace {
template <typename vertex_t, typename weight_t>
auto weaklyConnectedComponents(vertex_t numNodes,
                               boost::beast::span<const vertex_t>& rowIndices,
                               boost::beast::span<const vertex_t>& colIndices,
                               boost::beast::span<const weight_t>& edgeWeights,
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
       Acts::zip(rowIndices, colIndices, edgeWeights)) {
    boost::add_edge(row, col, weight, g);
  }

  return boost::connected_components(g, &trackLabels[0]);
}
}  // namespace

namespace Acts {

std::vector<std::vector<int>> BoostTrackBuilding::operator()(
    PipelineTensors tensors, std::vector<int>& spacepointIDs,
    const ExecutionContext& execContext) {
  ACTS_DEBUG("Start track building");

  using RTI = const Tensor<std::int64_t>&;
  const auto& edgeTensor =
      tensors.edgeIndex.device().isCpu()
          ? static_cast<RTI>(tensors.edgeIndex)
          : static_cast<RTI>(tensors.edgeIndex.clone(
                {Acts::Device::Cpu(), execContext.stream}));

  assert(tensors.scoreTensor.has_value());
  using RTF = const Tensor<float>&;
  const auto& scoreTensor =
      tensors.edgeScores->device().isCpu()
          ? static_cast<RTF>(*tensors.edgeScores)
          : static_cast<RTF>(tensors.edgeScores->clone(
                {Acts::Device::Cpu(), execContext.stream}));

  assert(edgeTensor.shape().at(0) == 2);
  assert(edgeTensor.shape().at(1) == scoreTensor.shape().at(0));

  const auto numSpacepoints = spacepointIDs.size();
  const auto numEdges = edgeTensor.shape().at(1);

  if (numEdges == 0) {
    ACTS_WARNING("No edges remained after edge classification");
    return {};
  }

  using vertex_t = std::int64_t;
  using weight_t = float;

  boost::beast::span<const vertex_t> rowIndices(edgeTensor.data(), numEdges);
  boost::beast::span<const vertex_t> colIndices(edgeTensor.data() + numEdges,
                                                numEdges);
  boost::beast::span<const weight_t> edgeWeights(scoreTensor.data(), numEdges);

  std::vector<vertex_t> trackLabels(numSpacepoints);

  auto numberLabels = weaklyConnectedComponents<vertex_t, weight_t>(
      numSpacepoints, rowIndices, colIndices, edgeWeights, trackLabels);

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
