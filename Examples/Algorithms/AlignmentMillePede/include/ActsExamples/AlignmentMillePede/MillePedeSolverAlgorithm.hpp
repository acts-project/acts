// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"

namespace ActsExamples {

/// @brief Algorithm for running a MillePede alignment fit at the end of the event loop.
/// This wraps a call to the MillePedeSolver.
class MillePedeSolverAlgorithm final : public IAlgorithm {
 public:
  /// configuration
  struct Config {
    ActsPlugins::MillePedeSolver::Config solverConfig;
  };

  /// Constructor of the sandbox algorithm
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  explicit MillePedeSolverAlgorithm(
      Config cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "MillePedeSolverAlgorithm", Acts::Logging::INFO))
      : IAlgorithm("MillePedeSolverAlgorithm", std::move(logger)),
        m_cfg(std::move(cfg)) {}

  /// Framework execute method of the sandbox algorithm
  ///
  /// For this algorithm, no work is done at event-loop time.
  /// Instead, the alignment fit will be executed in the finalize method.
  /// This allows to collect tracks over many events.
  ProcessCode execute(const AlgorithmContext& /*ctx*/) const override {
    return ProcessCode::SUCCESS;
  };
  /// finalize() method run after the event loop.
  /// Will call the Millepede solver and write a
  /// result file that can be picked up in subsequent
  /// runs.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// configuration instance
  Config m_cfg;
};

}  // namespace ActsExamples
