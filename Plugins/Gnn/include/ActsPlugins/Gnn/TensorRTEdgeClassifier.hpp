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
#include <vector>

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// Edge classifier using TensorRT inference
class TensorRTEdgeClassifier final : public EdgeClassificationBase {
 public:
  /// Configuration struct for TensorRT edge classifier
  struct Config {
    /// Path to the TensorRT model file
    std::string modelPath;
    /// List of feature indices to use for edge classification
    std::vector<int> selectedFeatures;
    /// Classification score threshold for edge filtering
    float cut = 0.5;

    /// Number of parallel execution contexts for inference
    std::size_t numExecutionContexts = 1;
  };

  /// Constructor
  /// @param cfg Configuration parameters
  /// @param logger Logging instance
  TensorRTEdgeClassifier(const Config &cfg,
                         std::unique_ptr<const Acts::Logger> logger);
  ~TensorRTEdgeClassifier();

  PipelineTensors operator()(PipelineTensors tensors,
                             const ExecutionContext &execContext = {}) override;

  /// Get the configuration
  /// @return Reference to the configuration struct
  const Config &config() const;

 private:
  class Impl;
  std::unique_ptr<Impl> m_impl;
};

/// @}
}  // namespace ActsPlugins
