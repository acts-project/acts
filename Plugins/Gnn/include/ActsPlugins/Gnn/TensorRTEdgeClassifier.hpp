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
#include <mutex>
#include <semaphore>
#include <vector>

namespace nvinfer1 {
class IRuntime;
class ICudaEngine;
class ILogger;
class IExecutionContext;
}  // namespace nvinfer1

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
  /// @return Copy of the configuration struct
  Config config() const { return m_cfg; }

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<nvinfer1::IRuntime> m_runtime;
  std::unique_ptr<nvinfer1::ICudaEngine> m_engine;
  std::unique_ptr<nvinfer1::ILogger> m_trtLogger;

  std::size_t m_maxEdges = 0;
  std::size_t m_maxNodes = 0;

  mutable std::optional<std::counting_semaphore<>> m_count;
  mutable std::mutex m_contextMutex;
  mutable std::vector<std::unique_ptr<nvinfer1::IExecutionContext>> m_contexts;
};

/// @}
}  // namespace ActsPlugins
