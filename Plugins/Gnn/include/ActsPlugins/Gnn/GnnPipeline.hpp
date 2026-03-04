// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Gnn/Stages.hpp"

#include <chrono>
#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Timing information for GNN pipeline execution
struct GnnTiming {
  /// Duration type in milliseconds
  using Duration = std::chrono::duration<float, std::milli>;

  /// Time spent building the graph
  Duration graphBuildingTime = Duration{0.f};
  /// Time spent in each classifier stage
  boost::container::small_vector<Duration, 3> classifierTimes;
  /// Time spent building tracks
  Duration trackBuildingTime = Duration{0.f};
};

/// Hook interface for GNN pipeline callbacks
class GnnHook {
 public:
  virtual ~GnnHook() = default;
  /// @brief Callback operator invoked during pipeline execution
  virtual void operator()(const PipelineTensors & /*tensors*/,
                          const ExecutionContext & /*execCtx*/) const {};
};

/// Graph Neural Network pipeline for track finding
class GnnPipeline {
 public:
  /// @brief Construct GNN pipeline
  /// @param graphConstructor Graph construction stage
  /// @param edgeClassifiers Edge classification stages
  /// @param trackBuilder Track building stage
  /// @param logger Logger instance
  GnnPipeline(
      std::shared_ptr<GraphConstructionBase> graphConstructor,
      std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
      std::shared_ptr<TrackBuildingBase> trackBuilder,
      std::unique_ptr<const Acts::Logger> logger);

  /// @brief Run the GNN pipeline
  /// @param features Input feature vector
  /// @param moduleIds Module identifiers for each space point
  /// @param spacePointIDs Space point identifiers
  /// @param device Device to run on (CPU/GPU)
  /// @param hook Optional callback hook for pipeline execution
  /// @param timing Optional timing output
  /// @return Vector of track candidates (each track is a vector of space point IDs)
  std::vector<std::vector<int>> run(std::vector<float> &features,
                                    const std::vector<std::uint64_t> &moduleIds,
                                    std::vector<int> &spacePointIDs,
                                    Device device, const GnnHook &hook = {},
                                    GnnTiming *timing = nullptr) const;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;

  std::shared_ptr<GraphConstructionBase> m_graphConstructor;
  std::vector<std::shared_ptr<EdgeClassificationBase>> m_edgeClassifiers;
  std::shared_ptr<TrackBuildingBase> m_trackBuilder;

  const Acts::Logger &logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
