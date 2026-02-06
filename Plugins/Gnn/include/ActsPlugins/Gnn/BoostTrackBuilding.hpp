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

/// Track building implementation using Boost
class BoostTrackBuilding final : public TrackBuildingBase {
 public:
  struct Config {};

  /// Constructor
  /// @param cfg Configuration object
  /// @param logger Logger instance
  BoostTrackBuilding(const Config &cfg,
                     std::unique_ptr<const Acts::Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  std::vector<std::vector<int>> operator()(
      PipelineTensors tensors, std::vector<int> &spacepointIDs,
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
