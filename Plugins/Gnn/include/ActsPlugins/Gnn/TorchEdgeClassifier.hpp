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

class TorchEdgeClassifier final : public EdgeClassificationBase {
 public:
  struct Config {
    std::string modelPath;
    std::vector<int> selectedFeatures = {};
    float cut = 0.5;
    int nChunks = 1;  // NOTE for GNN use 1
    bool undirected = false;
    int deviceID = 0;
    bool useEdgeFeatures = false;
  };

  TorchEdgeClassifier(const Config &cfg,
                      std::unique_ptr<const Acts::Logger> logger);
  ~TorchEdgeClassifier();

  PipelineTensors operator()(PipelineTensors tensors,
                             const ExecutionContext &execContext = {}) override;

  Config config() const { return m_cfg; }

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<torch::jit::Module> m_model;
};

/// @}
}  // namespace ActsPlugins
