// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"
#include "ActsPlugins/Mille/MillePedeSteering.hpp"

namespace ActsExamples {

/// @brief Algorithm for running a MillePede alignment fit at the end of the event loop.
/// This wraps a call to the MillePedeSolver.
class MillePedeSolverAlgorithm final : public IAlgorithm {
 public:
  /// configuration
  struct Config {
    ActsPlugins::ActsToMille::MillePedeSolver::Config solverConfig;
    ActsPlugins::ActsToMille::MillePedeSteering::Config steeringConfig;
  };

  /// Constructor of the sandbox algorithm
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  explicit MillePedeSolverAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr)
      : IAlgorithm("ActsSolverFromMille", std::move(logger)),
        m_cfg(std::move(cfg)) {}

  /// Framework execute method of the sandbox algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  /// Nothing to do here
  ProcessCode execute(const AlgorithmContext& /*ctx*/) const override {
    return ProcessCode::SUCCESS;
  };
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// configuration instance
  Config m_cfg;
};

}  // namespace ActsExamples
