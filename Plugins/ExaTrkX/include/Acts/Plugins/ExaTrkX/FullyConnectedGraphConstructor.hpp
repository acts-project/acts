// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// Graph construction class that creates a fully connected graph
/// Mostly for debugging purposes as graph size is O(n^2)
class FullyConnectedGraphConstructor : public GraphConstructionBase {
 public:
  struct Config {
    // Due to O(n^2) the algorithm can be configured to fail above a certain
    // graph size.
    std::size_t maxGraphSize = 1000;

    // Only connect edges up to this delta R
    float maxDeltaR = std::numeric_limits<float>::infinity();

    // Offset of the r feature in the input features
    std::size_t rOffset = 0;

    // Scales (used for edge feature building)
    float rScale = 1.0;
    float phiScale = 1.0;
    float zScale = 1.0;
    float etaScale = 1.0;
  };

  FullyConnectedGraphConstructor(const Config &cfg,
                                 std::unique_ptr<const Acts::Logger> _logger)
      : m_cfg(cfg), m_logger(std::move(_logger)) {}

  std::tuple<std::any, std::any, std::any> operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<std::uint64_t> &moduleIds,
      const ExecutionContext &execContext = {}) override;

  torch::Device device() const override {
    return torch::Device(torch::kCUDA, 0);
  }

  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  const auto &logger() const { return *m_logger; }
};

}  // namespace Acts
