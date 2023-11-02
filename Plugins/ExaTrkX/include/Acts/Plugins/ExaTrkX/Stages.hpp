// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <vector>

namespace Acts {

// TODO maybe replace std::any with some kind of variant<unique_ptr<torch>,
// unique_ptr<onnx>>?
// TODO maybe replace input for GraphConstructionBase with some kind of
// boost::multi_array / Eigen::Array

class GraphConstructionBase {
 public:
  /// Perform the graph construction
  ///
  /// @param inputValues Flattened input data
  /// @param numNodes Number of nodes. inputValues.size() / numNodes
  /// then gives the number of features
  /// @param deviceHint Which GPU to pick. Not relevant for CPU-only builds
  ///
  /// @return (node_tensor, edge_tensore)
  virtual std::tuple<std::any, std::any> operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      int deviceHint = -1) = 0;

  virtual ~GraphConstructionBase() = default;
};

class EdgeClassificationBase {
 public:
  /// Perform edge classification
  ///
  /// @param nodes Node tensor with shape (n_nodes, n_node_features)
  /// @param edges Edge-index tensor with shape (2, n_edges)
  /// @param deviceHint Which GPU to pick. Not relevant for CPU-only builds
  ///
  /// @return (node_tensor, edge_tensor, score_tensor)
  virtual std::tuple<std::any, std::any, std::any> operator()(
      std::any nodes, std::any edges, int deviceHint = -1) = 0;

  virtual ~EdgeClassificationBase() = default;
};

class TrackBuildingBase {
 public:
  /// Perform track building
  ///
  /// @param nodes Node tensor with shape (n_nodes, n_node_features)
  /// @param edges Edge-index tensor with shape (2, n_edges)
  /// @param edgeWeights Edge-weights of the previous edge classification phase
  /// @param spacepointIDs IDs of the nodes (must have size=n_nodes)
  /// @param deviceHint Which GPU to pick. Not relevant for CPU-only builds
  ///
  /// @return tracks (as vectors of node-IDs)
  virtual std::vector<std::vector<int>> operator()(
      std::any nodes, std::any edges, std::any edgeWeights,
      std::vector<int> &spacepointIDs, int deviceHint = -1) = 0;

  virtual ~TrackBuildingBase() = default;
};

}  // namespace Acts
