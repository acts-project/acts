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

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4VUserDetectorConstruction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4MagneticField;
class G4VUserPhysicsList;
class G4FieldManager;

namespace Acts {
class TrackingGeometry;
class MagneticFieldProvider;
class Volume;
}  // namespace Acts

namespace ActsExamples {

class DetectorConstructionFactory;
class SensitiveSurfaceMapper;
class EventStoreHolder;
struct Geant4Instance;

class Geant4SimulationBase : public IAlgorithm {
 public:
  /// Nested configuration struct for the Geant4 simulation
  struct Config {
    // Name of the input particle collection
    std::string inputParticles = "";

    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;

    /// Detector construction object.
    /// G4RunManager will take care of deletion
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory;
  };

  Geant4SimulationBase(const Config& cfg, std::string name,
                       Acts::Logging::Level level = Acts::Logging::INFO);

  ~Geant4SimulationBase() override;

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const override;

  /// Readonly access to the configuration
  virtual const Config& config() const = 0;

 protected:
  void initializeCommon(const Config& cfg, G4VUserPhysicsList* physicsList);
  void kickRunManager();

  std::unique_ptr<const Acts::Logger> m_logger;

  std::shared_ptr<EventStoreHolder> m_eventStoreHolder;

  std::shared_ptr<Geant4Instance> m_gean4Instance;

  /// Our Geant4Manager is taking care of the lifetime
  G4RunManager* m_runManager{};

  /// The G4 physics list
  G4VUserPhysicsList* m_physicsList{};

  /// Detector construction object.
  /// G4RunManager will take care of deletion
  G4VUserDetectorConstruction* m_detectorConstruction{};

  /// User Action: Primary generator action of the simulation
  /// G4RunManager will take care of deletion
  G4VUserPrimaryGeneratorAction* m_primaryGeneratorAction{};

  /// User Action: Run
  /// G4RunManager will take care of deletion
  G4UserRunAction* m_runAction{};

  /// User Action: Event
  /// G4RunManager will take care of deletion
  G4UserEventAction* m_eventAction{};

  /// User Action: Tracking
  /// G4RunManager will take care of deletion
  G4UserTrackingAction* m_trackingAction{};

  /// User Action: Stepping
  /// G4RunManager will take care of deletion
  G4UserSteppingAction* m_steppingAction{};

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
};

/// Algorithm to run Geant4 simulation in the ActsExamples framework
///
/// This algorithm can be configured with the standardard Geant4
/// components:
/// - (a) the generator Action
/// - (b) detector construction and magnetic field
/// - (c) the user actions
///
/// In order to run within the ACTS framework, acces to the
/// EventData is provided by a EventStoreHolder which provides
/// individual slots for the event containers and the store.
///
/// The Geant4Simulation algorithm clears those after processing.
class Geant4Simulation final : public Geant4SimulationBase {
 public:
  struct Config : public Geant4SimulationBase::Config {
    // Name of the output collection : hits
    std::string outputSimHits = "simhits";

    // Name of the output collection : initial particles
    std::string outputParticlesInitial = "particles_initial";

    // Name of the output collection : final particles
    std::string outputParticlesFinal = "particles_final";

    /// The ACTS tracking geometry
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// The ACTS Magnetic field provider
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    std::vector<std::string> volumeMappings = {"Silicon"};

    std::vector<std::string> materialMappings;

    std::shared_ptr<const Acts::Volume> killVolume;

    double killAfterTime = std::numeric_limits<double>::infinity();

    bool recordHitsOfSecondaries = true;

    bool keepParticlesWithoutHits = true;
  };

  /// Simulation constructor
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
  const Config& config() const final { return m_cfg; }

 private:
  Config m_cfg;

  /// The (wrapped) ACTS Magnetic field provider as a Geant4 module
  std::unique_ptr<G4MagneticField> m_magneticField;
  std::unique_ptr<G4FieldManager> m_fieldManager;

  WriteDataHandle<SimParticleContainer> m_outputParticlesInitial{
      this, "OutputParticlesInitial"};
  WriteDataHandle<SimParticleContainer> m_outputParticlesFinal{
      this, "OutputParticlesFinal"};
  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHIts"};
};

class Geant4MaterialRecording final : public Geant4SimulationBase {
 public:
  struct Config : public Geant4SimulationBase::Config {
    // Name of the output collection: material tracks
    std::string outputMaterialTracks = "material_tracks";
  };

  /// Material recording constructor
  ///
  /// @param config is the configuration struct
  /// @param level is the logging level to be used
  Geant4MaterialRecording(const Config& config,
                          Acts::Logging::Level level = Acts::Logging::INFO);

  ~Geant4MaterialRecording() override;

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  /// Readonly access to the configuration
  const Config& config() const final { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<std::unordered_map<size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};
};

}  // namespace ActsExamples
