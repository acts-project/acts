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

#include <memory>

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Track building stage based on the EdgeLayerConnector algorithm from
/// the ModuleMapGraph library (CUDA based)
class EdgeLayerConnector final : public TrackBuildingBase {
 public:
  /// Configuration for the EdgeLayerConnector
  struct Config {
    /// Number of thread blocks for parallel edge processing
    std::size_t blockSize = 512;
    /// Maximum number of hits allowed per track candidate
    std::size_t maxHitsPerTrack = 30;
    /// Minimum number of hits required to keep a track candidate
    std::size_t minHits = 3;
    /// Edge weight threshold below which edges are discarded
    float weightsCut = 0.01;
  };

  /// @param cfg Configuration struct
  /// @param logger Logger instance
  EdgeLayerConnector(const Config &cfg,
                     std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  std::vector<std::vector<int>> operator()(
      PipelineTensors tensors, std::vector<int> &spacepointIDs,
      const ExecutionContext &execContext = {}) override;

  /// @return Read-only reference to the configuration
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
