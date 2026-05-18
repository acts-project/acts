// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/Tensor.hpp"

#include <cstdint>
#include <exception>
#include <optional>
#include <vector>

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Error that is thrown if no edges are found
struct NoEdgesError : std::exception {};

/// Struct that ties together the tensors used in the GNN pipeline
struct PipelineTensors {
  /// Tensor containing node feature data
  Tensor<float> nodeFeatures;
  /// Tensor containing edge connectivity indices
  Tensor<std::int64_t> edgeIndex;
  /// Optional tensor containing edge feature data
  std::optional<Tensor<float>> edgeFeatures;
  /// Optional tensor containing edge classification scores
  std::optional<Tensor<float>> edgeScores;
};

/// Base class for graph construction algorithms
class GraphConstructionBase {
 public:
  /// Perform the graph construction
  ///
  /// @param inputValues Flattened input data
  /// @param numNodes Number of nodes. inputValues.size() / numNodes
  /// then gives the number of features
  /// @param moduleIds Module IDs of the features (used for module-map-like
  /// graph construction)
  /// @param execContext Device & stream information
  /// @return (node_features, edge_features, edge_index)
  virtual PipelineTensors operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<std::uint64_t> &moduleIds,
      const ExecutionContext &execContext = {}) = 0;

  virtual ~GraphConstructionBase() = default;
};

/// @brief Base class for edge classification stage in GNN pipeline
class EdgeClassificationBase {
 public:
  /// Perform edge classification
  ///
  /// @param tensors Input pipeline tensors
  /// @param execContext Device & stream information
  ///
  /// @return (node_features, edge_features, edge_index, edge_scores)
  virtual PipelineTensors operator()(
      PipelineTensors tensors, const ExecutionContext &execContext = {}) = 0;

  virtual ~EdgeClassificationBase() = default;
};

/// Base class for track building implementations
class TrackBuildingBase {
 public:
  /// Perform track building
  ///
  /// @param tensors Input pipeline tensors
  /// @param spacePointIDs IDs of the nodes (must have size=n_nodes)
  /// @param execContext Device & stream information
  ///
  /// @return tracks (as vectors of node-IDs)
  virtual std::vector<std::vector<int>> operator()(
      PipelineTensors tensors, std::vector<int> &spacePointIDs,
      const ExecutionContext &execContext = {}) = 0;

  virtual ~TrackBuildingBase() = default;
};

/// @}
}  // namespace ActsPlugins
