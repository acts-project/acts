// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <memory>
#include <mutex>
#include <string>

#include "G4VUserDetectorConstruction.hh"

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4MagneticField;
class G4VUserPhysicsList;

namespace ActsExamples {

class SensitiveSurfaceMapper;

/// Algorithm to run Geant4 simulation in the ActsExamples framework
///
/// This algorithm can be configured with the standardard Geant4
/// components:
/// - (a) the generator Action
/// - (b) detector construction and magnetic field
/// - (c) the user actions
///
/// In order to run within the ACTS framework, acces to the
/// EventData is provided by a EventStoreRegistry which provides
/// individual slots for the event containers and the store.
///
/// The Geant4Simulation algorithm clears those after processing.
class Geant4Simulation final : public IAlgorithm {
 public:
  /// Nested configuration struct for the Geant4 simulation
  struct Config {
    // Name of the input particle collection
    std::string inputParticles = "";

    // Name of the output collection : hits
    std::string outputSimHits = "";

    // Name of the output collection : initial particles
    std::string outputParticlesInitial = "";

    // Name of the output collection : final particles
    std::string outputParticlesFinal = "";

    // Name of the output collection: material tracks
    std::string outputMaterialTracks = "";

    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;

    /// The G4 run manager
    std::shared_ptr<G4RunManager> runManager;

    /// User Action: Primary generator action of the simulation
    G4VUserPrimaryGeneratorAction* primaryGeneratorAction = nullptr;

    /// User Action: Run
    G4UserRunAction* runAction = nullptr;

    /// User Action: Event
    G4UserEventAction* eventAction = nullptr;

    /// User Action: Tracking
    G4UserTrackingAction* trackingAction = nullptr;

    /// User Action: Stepping
    G4UserSteppingAction* steppingAction = nullptr;

    /// Detector construction object.
    G4VUserDetectorConstruction* detectorConstruction = nullptr;

    /// The (wrapped) ACTS Magnetic field provider as a Geant4 module
    G4MagneticField* magneticField = nullptr;

    // The ACTS to Geant4 sensitive wrapper
    std::shared_ptr<const SensitiveSurfaceMapper> sensitiveSurfaceMapper =
        nullptr;
  };

  /// Constructor with arguments
  ///
  /// @param config is the configuration struct
  /// @param level is the logging level to be used
  Geant4Simulation(const Config& config,
                   Acts::Logging::Level level = Acts::Logging::INFO);
  ~Geant4Simulation() override;

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  /// Readonly access to the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  WriteDataHandle<SimParticleContainer> m_outputParticlesInitial{
      this, "OutputParticlesInitial"};
  WriteDataHandle<SimParticleContainer> m_outputParticlesFinal{
      this, "OutputParticlesFinal"};
  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHIts"};
  WriteDataHandle<std::unordered_map<size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};

  // Has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};

}  // namespace ActsExamples
