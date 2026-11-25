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
#include <functional>
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
  /// @param shrinkNodes If true, unused nodes are removed from the graph
  /// before each edge classification stage
  GnnPipeline(
      std::shared_ptr<GraphConstructionBase> graphConstructor,
      std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
      std::shared_ptr<TrackBuildingBase> trackBuilder,
      std::unique_ptr<const Acts::Logger> logger, bool shrinkNodes = false);

  /// @brief Construct GNN pipeline with a seeded graph construction stage
  /// @param seededGraphConstructor Graph construction stage that takes seeds
  /// @param edgeClassifiers Edge classification stages
  /// @param trackBuilder Track building stage
  /// @param logger Logger instance
  /// @param shrinkNodes If true, unused nodes are removed from the graph
  /// before each edge classification stage
  GnnPipeline(
      std::shared_ptr<SeededGraphConstructionBase> seededGraphConstructor,
      std::vector<std::shared_ptr<EdgeClassificationBase>> edgeClassifiers,
      std::shared_ptr<TrackBuildingBase> trackBuilder,
      std::unique_ptr<const Acts::Logger> logger, bool shrinkNodes = false);

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

  /// @brief Run the GNN pipeline with seeds
  /// Requires a pipeline built with a seeded graph construction stage.
  /// @param features Input feature vector
  /// @param seeds Seeds passed to the graph construction stage. Seed hit
  /// indices refer to the nodes of @p features
  /// @param moduleIds Module identifiers for each space point
  /// @param spacePointIDs Space point identifiers
  /// @param device Device to run on (CPU/GPU)
  /// @param hook Optional callback hook for pipeline execution
  /// @param timing Optional timing output
  /// @return Vector of track candidates (each track is a vector of space point IDs)
  std::vector<std::vector<int>> run(std::vector<float> &features,
                                    const std::vector<Seed> &seeds,
                                    const std::vector<std::uint64_t> &moduleIds,
                                    std::vector<int> &spacePointIDs,
                                    Device device, const GnnHook &hook = {},
                                    GnnTiming *timing = nullptr) const;

 private:
  using GraphBuilder = std::function<PipelineTensors(const ExecutionContext &)>;

  /// Shared pipeline body; buildGraph runs the graph construction stage
  std::vector<std::vector<int>> runPipeline(const GraphBuilder &buildGraph,
                                            std::vector<int> &spacePointIds,
                                            Device device, const GnnHook &hook,
                                            GnnTiming *timing) const;

  /// Throws if the stages are not consistently configured
  void validateStages() const;

  std::unique_ptr<const Acts::Logger> m_logger;

  std::shared_ptr<GraphConstructionBase> m_graphConstructor;
  std::shared_ptr<SeededGraphConstructionBase> m_seededGraphConstructor;
  std::vector<std::shared_ptr<EdgeClassificationBase>> m_edgeClassifiers;
  std::shared_ptr<TrackBuildingBase> m_trackBuilder;
  bool m_shrinkNodes = false;

  const Acts::Logger &logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
