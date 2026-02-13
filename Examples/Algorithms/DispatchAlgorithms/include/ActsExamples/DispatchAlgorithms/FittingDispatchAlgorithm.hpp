// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

namespace ActsExamples {

/// This Algorithm retrieves Acts event data for pattern recognition 
/// or re-fitting and forwards it to a user defined function for further processing
class PatternDispatchAlgorithm final : public IAlgorithm {
 public:

 /// Configuration class it allows to connect to python functions
 class Config {
   public:
   

  };

  /// Construct the smearing algorithm.
  ///
  /// @param config is the algorithm configuration
  /// @param level is the logging level
  PatternDispatchAlgorithm(Config config, Acts::Logging::Level level);

  /// Build measurement from simulation hits at input.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 private:

  /// Configuration of the Algorithm
  Config m_cfg;


};

}  // namespace ActsExamples
