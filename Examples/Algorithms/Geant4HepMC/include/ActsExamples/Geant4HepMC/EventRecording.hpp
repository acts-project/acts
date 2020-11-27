// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <mutex>
#include <G4VUserDetectorConstruction.hh>
#include <HepMC3/GenEvent.h>

class G4RunManager;

namespace ActsExamples {

class EventRecording final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    /// The input collection of particles
    std::string inputParticles = "";
    /// The recorded events output
    std::string outputHepMcTracks = "geant-outcome-tracks";

    std::unique_ptr<G4VUserDetectorConstruction> detectorConstruction = nullptr;

    /// random number seed 1
    int seed1 = 12345;
    /// random number seed 2
    int seed2 = 45678;

    /// List of processes that can be combined to a single vertex
    std::vector<std::string> processesCombine;
    /// Optional selective recording based on a process
    /// @note All events are recorded if this is empty
    std::string processSelect;
    /// List to veto events with certain processes
    std::vector<std::string> processesReject;
  };

  /// Constructor
  EventRecording(Config&& cnf, Acts::Logging::Level level);
  ~EventRecording();

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
  /// G4 run manager
  std::unique_ptr<G4RunManager> m_runManager;

  // has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};
}  // namespace ActsExamples
