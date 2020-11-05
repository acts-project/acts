// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

class G4RunManager;

namespace ActsExamples {

struct EventFraction {
  EventFraction() = default;

  EventFraction(ActsExamples::SimParticle& initPart,
                std::vector<ActsExamples::SimParticle>& finalPart)
      : initialParticle(initPart), finalParticles(std::move(finalPart)) {}

  ActsExamples::SimParticle initialParticle;
  std::vector<ActsExamples::SimParticle> finalParticles;

  bool soft = false;
  unsigned int multiplicity = 0;
};

/// @brief This class extracts a certain process from a HepMC event record.
class EventExtraction final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
	  /// The input collection
    std::string inputEvents;
    /// The output collection
    std::string outputEventFraction = "event-fraction";
    /// The process that should be extracted
    std::string extractionProcess;
    
    /// Minimum absolute value of considered PDG IDs
    int minAbsPdg = 40;
    /// Maximum absolute value of considered PDG IDs
    int maxAbsPdg = 2212;
    /// Minimum momentum of considered particles
    double pMin = 50. * Acts::UnitConstants::MeV;
  };

  /// Constructor
  EventExtraction(Config&& cnf, Acts::Logging::Level level);
  ~EventExtraction();

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
};
}  // namespace ActsExamples
