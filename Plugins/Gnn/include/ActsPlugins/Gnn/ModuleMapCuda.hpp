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
#include <string>

namespace ActsPlugins {

/// @addtogroup gnn_plugin
/// @{

/// CUDA-based module map for graph construction
class ModuleMapCuda : public GraphConstructionBase {
 public:
  /// Configuration for ModuleMapCuda
  struct Config {
    /// Path to module map file
    std::string moduleMapPath;
    /// Radial coordinate scaling factor
    float rScale = 1.0;
    /// Azimuthal coordinate scaling factor
    float phiScale = 1.0;
    /// Z-coordinate scaling factor
    float zScale = 1.0;
    /// Pseudorapidity scaling factor
    float etaScale = 1.0;

    /// Enable more parallel execution
    bool moreParallel = true;
    /// CUDA device ID
    int gpuDevice = 0;
    /// Number of GPU blocks
    int gpuBlocks = 512;

    /// Small numerical constant for stability
    float epsilon = 1e-8f;
  };

  /// Constructor
  /// @param cfg Configuration parameters
  /// @param logger Logger instance
  ModuleMapCuda(const Config &cfg, std::unique_ptr<const Acts::Logger> logger);

  ~ModuleMapCuda() override;

  /// Access configuration
  /// @return Configuration reference
  const auto &config() const { return m_cfg; }

  PipelineTensors operator()(std::vector<float> &inputValues,
                             std::size_t numNodes,
                             const std::vector<std::uint64_t> &moduleIds,
                             const ExecutionContext &execContext = {}) override;

 private:
  class Impl;
  std::unique_ptr<Impl> m_impl;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  const auto &logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
