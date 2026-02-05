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

namespace Ort {
class Env;
class Session;
class Value;
}  // namespace Ort

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

/// ONNX-based edge classifier for GNN pipeline
class OnnxEdgeClassifier final : public EdgeClassificationBase {
 public:
  /// Configuration for OnnxEdgeClassifier
  struct Config {
    /// Path to the ONNX model file
    std::string modelPath;
    /// Classification threshold cut
    float cut = 0.5;
  };

  /// Constructor
  /// @param cfg Configuration parameters
  /// @param logger Logger instance
  OnnxEdgeClassifier(const Config &cfg,
                     std::unique_ptr<const Acts::Logger> logger);
  ~OnnxEdgeClassifier();

  PipelineTensors operator()(PipelineTensors tensors,
                             const ExecutionContext &execContext = {}) override;

  /// Get configuration
  /// @return Configuration parameters
  Config config() const { return m_cfg; }

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<Ort::Env> m_env;
  std::unique_ptr<Ort::Session> m_model;

  std::vector<std::string> m_inputNames;
  std::string m_outputName;
};

/// @}
}  // namespace ActsPlugins
