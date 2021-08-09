// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ExtractedSimulationProcess.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

class G4RunManager;

namespace ActsExamples {

/// @brief This class extracts a certain process from a HepMC event record.
class HepMCProcessExtractor final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    /// The input collection
    std::string inputEvents;
    /// The output collection
    std::string outputSimulationProcesses = "event-fraction";
    /// The process that should be extracted
    std::string extractionProcess;

    /// Minimum absolute value of considered PDG IDs
    int absPdgMin = 40;
    /// Maximum absolute value of considered PDG IDs
    int absPdgMax = 2212;
    /// Minimum momentum of considered particles
    double pMin = 50. * Acts::UnitConstants::MeV;
  };

  /// Constructor
  /// @param config the configuration
  /// @param level the log level
  HepMCProcessExtractor(Config config, Acts::Logging::Level level);
  ~HepMCProcessExtractor();

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
};

}  // namespace ActsExamples
