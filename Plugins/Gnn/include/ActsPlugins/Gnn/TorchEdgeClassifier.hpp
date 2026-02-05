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

/// Edge classifier using PyTorch inference
class TorchEdgeClassifier final : public EdgeClassificationBase {
 public:
  /// Configuration struct for Torch edge classifier
  struct Config {
    /// Path to the PyTorch model file
    std::string modelPath;
    /// Selected feature indices for input
    std::vector<int> selectedFeatures = {};
    /// Classification score threshold for edge filtering
    float cut = 0.5;
    /// Number of chunks to process
    int nChunks = 1;  // NOTE for GNN use 1
    /// Whether to treat graph as undirected
    bool undirected = false;
    /// CUDA device ID to use for inference
    int deviceID = 0;
    /// Whether to use edge features
    bool useEdgeFeatures = false;
  };

  /// Constructor
  /// @param cfg Configuration parameters
  /// @param logger Logging instance
  TorchEdgeClassifier(const Config &cfg,
                      std::unique_ptr<const Acts::Logger> logger);
  ~TorchEdgeClassifier();

  PipelineTensors operator()(PipelineTensors tensors,
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
