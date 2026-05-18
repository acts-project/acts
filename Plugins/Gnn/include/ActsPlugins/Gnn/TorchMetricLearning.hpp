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

/// @cond
namespace torch::jit {
class Module;
}

namespace c10 {
enum class DeviceType : std::int8_t;
}
/// @endcond

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Graph construction using PyTorch metric learning
class TorchMetricLearning final : public GraphConstructionBase {
 public:
  /// Configuration struct for Torch metric learning
  struct Config {
    /// Path to the PyTorch model file
    std::string modelPath;
    /// Selected feature indices for input
    std::vector<int> selectedFeatures = {};
    /// Dimensionality of the embedding space
    int embeddingDim = 8;
    /// Radius value for graph construction
    float rVal = 1.6;
    /// Number of nearest neighbors
    int knnVal = 500;
    /// Whether to shuffle edge directions
    bool shuffleDirections = false;
    /// CUDA device ID to use for inference
    int deviceID = 0;  // default is the first GPU if available

    /// Scaling factor for phi coordinate in edge features
    float phiScale = 3.141592654;
  };

  /// Constructor
  /// @param cfg Configuration parameters
  /// @param logger Logging instance
  TorchMetricLearning(const Config &cfg,
                      std::unique_ptr<const Acts::Logger> logger);
  ~TorchMetricLearning();

  PipelineTensors operator()(std::vector<float> &inputValues,
                             std::size_t numNodes,
                             const std::vector<std::uint64_t> &moduleIds,
                             const ExecutionContext &execContext = {}) override;

  /// Get the configuration
  /// @return Copy of the configuration struct
  Config config() const { return m_cfg; }

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<torch::jit::Module> m_model;
};

/// @}
}  // namespace ActsPlugins
