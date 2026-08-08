// This file is part of the ACTS project
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Gnn/Stages.hpp"

#include <memory>
#include <string>

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Track building stage using the D-WALK graph segmentation algorithm.
///
/// The initial connected-component score cut is expected to be applied by the
/// edge-classifier stage. This stage consumes the already-filtered graph and
/// applies the D-WALK simple-path and dynamic-programming steps.
class DWalkTrackBuilding final : public TrackBuildingBase {
 public:
  /// Configuration for D-WALK track building
  struct Config {
    /// Minimal score for the best incoming/outgoing edge in the DP max-add cut
    float thMin = 0.1;
    /// Score above which all incoming/outgoing edges are kept in the DP max-add cut
    float thAdd = 0.6;
    /// Minimum number of hits required to keep a track candidate
    std::size_t minCandidateSize = 3;
    /// Node-feature column used for radial ordering
    std::size_t radialFeatureIndex = 0;
    /// DP path metric: "score_weighted_length" or "length"
    std::string pathMetric = "score_weighted_length";
  };

  /// Constructor
  /// @param cfg Configuration object
  /// @param logger Logger instance
  DWalkTrackBuilding(const Config &cfg,
                     std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  std::vector<std::vector<int>> operator()(
      PipelineTensors tensors, std::vector<int> &spacePointIDs,
      const ExecutionContext &execContext = {}) override;

  /// Get configuration
  /// @return Configuration object
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
