// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace torch::jit {
class Module;
}

namespace c10 {
enum class DeviceType : int8_t;
}

namespace Acts {

class TorchMetricLearning final : public Acts::GraphConstructionBase {
 public:
  struct Config {
    std::string modelPath;
    int numFeatures = 3;
    int embeddingDim = 8;
    float rVal = 1.6;
    int knnVal = 500;
    bool shuffleDirections = false;
  };

  TorchMetricLearning(const Config &cfg, std::unique_ptr<const Logger> logger);
  ~TorchMetricLearning();

  std::tuple<std::any, std::any> operator()(std::vector<float> &inputValues,
                                            std::size_t numNodes,
                                            int deviceHint = -1) override;

  Config config() const { return m_cfg; }

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;
  c10::DeviceType m_deviceType;
  std::unique_ptr<torch::jit::Module> m_model;
};

}  // namespace Acts
