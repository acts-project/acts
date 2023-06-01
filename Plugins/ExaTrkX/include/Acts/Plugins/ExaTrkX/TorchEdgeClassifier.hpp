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

class TorchEdgeClassifier final : public Acts::EdgeClassificationBase {
 public:
  struct Config {
    std::string modelPath;
    float cut = 0.21;
    int nChunks = 1;  // NOTE for GNN use 1
  };

  TorchEdgeClassifier(const Config &cfg);
  ~TorchEdgeClassifier();

  std::tuple<std::any, std::any, std::any> operator()(
      std::any nodes, std::any edges, const Logger &logger) override;

  Config config() const { return m_cfg; }

 private:
  Config m_cfg;
  c10::DeviceType m_deviceType;
  std::unique_ptr<torch::jit::Module> m_model;
};

}  // namespace Acts
