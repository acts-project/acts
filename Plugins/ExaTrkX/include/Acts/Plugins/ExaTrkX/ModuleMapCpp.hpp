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
#include <string>

template <typename T>
class graph_creator;

namespace Acts {

class ModuleMapCpp : public GraphConstructionBase {
 public:
  struct Config {
    std::string moduleMapPath;
  };

 private:
  Config m_cfg;
  std::unique_ptr<graph_creator<float>> m_graphCreator;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::vector<uint64_t> m_uniqueDoupletModuleIds;

  const auto &logger() const { return *m_logger; }

 public:
  ModuleMapCpp(const Config &cfg, std::unique_ptr<const Acts::Logger> logger);
  ~ModuleMapCpp();

  const auto &config() const { return m_cfg; }

  std::tuple<std::any, std::any, std::any> operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<uint64_t> &moduleIds,
      torch::Device device = torch::Device(torch::kCPU)) override;

  torch::Device device() const override { return torch::Device(torch::kCPU); }
};

}  // namespace Acts
