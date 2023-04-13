// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <iostream>
#include <stdexcept>

#include <G4EmParameters.hh>
#include <G4FieldManager.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4TransportationManager.hh>
#include <G4UniformMagField.hh>
#include <G4UserEventAction.hh>
#include <G4UserLimits.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4Version.hh>

namespace {
/// Helper method to add the user actions
/// @tparam manager_t the run manager type
/// @tparam actions_t the actions iterable list
template <typename manager_t, typename actions_t>
void setUserActions(manager_t& manager, actions_t& actions) {
  for (const auto& action : actions) {
    if (action != nullptr) {
      manager.SetUserAction(action);
    }
  }
}
}  // namespace

ActsExamples::Geant4Simulation::Geant4Simulation(
    const ActsExamples::Geant4Simulation::Config& config,
    Acts::Logging::Level level)
    : IAlgorithm("Geant4Simulation", level), m_cfg(config) {
  if (m_cfg.detectorConstruction == nullptr) {
    throw std::invalid_argument("Missing G4 DetectorConstruction object");
  }
  if (m_cfg.primaryGeneratorAction == nullptr) {
    throw std::invalid_argument("Missing G4 PrimaryGeneratorAction object");
  }
  if (!m_cfg.runManager) {
    throw std::invalid_argument("Missing G4 RunManager object");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  if (m_cfg.sensitiveSurfaceMapper) {
    if (m_cfg.outputSimHits.empty()) {
      ACTS_WARNING("No output sim hits collection configured");
    }
    m_outputSimHits.initialize(m_cfg.outputSimHits);

    if (m_cfg.outputParticlesInitial.empty()) {
      ACTS_WARNING("No output initial particles collection configured");
    }
    m_outputParticlesInitial.initialize(m_cfg.outputParticlesInitial);

    if (m_cfg.outputParticlesFinal.empty()) {
      ACTS_WARNING("No output final particles collection configured");
    }
    m_outputParticlesFinal.initialize(m_cfg.outputParticlesFinal);
  }

  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }

  // If we are in VERBOSE mode, set the verbose level in Geant4 to 2.
  // 3 would be also possible, but that produces infinite amount of output.
  const int geantVerboseLevel =
      logger().level() == Acts::Logging::VERBOSE ? 2 : 0;
  m_cfg.runManager->SetVerboseLevel(geantVerboseLevel);
  G4EventManager::GetEventManager()->SetVerboseLevel(geantVerboseLevel);
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(
      geantVerboseLevel);
  G4EventManager::GetEventManager()->GetStackManager()->SetVerboseLevel(
      geantVerboseLevel);

  // Suppress the printing of physics information.
#if G4VERSION_NUMBER >= 1100
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4EmParameters::Instance()->SetIsPrintedFlag(true);
#endif

  // Set the detector construction
  m_cfg.runManager->SetUserInitialization(m_cfg.detectorConstruction);

  // Set the primary generator action
  m_cfg.runManager->SetUserAction(m_cfg.primaryGeneratorAction);

  // Set the configured user actions
  setUserActions(*m_cfg.runManager, m_cfg.runActions);
  setUserActions(*m_cfg.runManager, m_cfg.eventActions);
  setUserActions(*m_cfg.runManager, m_cfg.trackingActions);
  setUserActions(*m_cfg.runManager, m_cfg.steppingActions);

  // Initialize the Geant4 run manager
  m_cfg.runManager->Initialize();

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Set the magnetic field
  if (m_cfg.magneticField != nullptr) {
    ACTS_INFO("Setting ACTS configured field to Geant4.");
    // Get the g4World cache
    G4VPhysicalVolume* g4World = m_cfg.detectorConstruction->Construct();
    /// Set the field ot the G4Field manager
    G4FieldManager* fieldMgr = new G4FieldManager();
    fieldMgr->SetDetectorField(m_cfg.magneticField);
    fieldMgr->CreateChordFinder(m_cfg.magneticField);

    // Propagate down to all childrend
    g4World->GetLogicalVolume()->SetFieldManager(fieldMgr, true);
  }

  // Map simulation to reconstruction geometry
  if (m_cfg.sensitiveSurfaceMapper != nullptr) {
    ACTS_INFO(
        "Remapping selected volumes from Geant4 to Acts::Surface::GeometryID");

    G4VPhysicalVolume* g4World = m_cfg.detectorConstruction->Construct();
    int sCounter = 0;
    m_cfg.sensitiveSurfaceMapper->remapSensitiveNames(
        g4World, Acts::Transform3::Identity(), sCounter);

    ACTS_INFO("Remapping successful for " << sCounter << " selected volumes.");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputMaterialTracks.maybeInitialize(m_cfg.outputMaterialTracks);
}

ActsExamples::Geant4Simulation::~Geant4Simulation() = default;

ActsExamples::ProcessCode ActsExamples::Geant4Simulation::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Ensure exclusive access to the Geant4 run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Set the seed new per event, so that we get reproducible results
  G4Random::setTheSeed(m_cfg.randomNumbers->generateSeed(ctx));

  // Get and reset event registry state
  auto& eventData = EventStoreRegistry::eventData();
  eventData = EventStoreRegistry::State{};

  // Register the current event store to the registry
  // this will allow access from the User*Actions
  eventData.store = &(ctx.eventStore);

  // Register the input particle read handle
  eventData.inputParticles = &m_inputParticles;

  ACTS_DEBUG("Sending Geant RunManager the BeamOn() command.");
  // Start simulation. each track is simulated as a separate Geant4 event.
  m_cfg.runManager->BeamOn(1);

  // Output handling: Initial/Final particles
  if (not m_cfg.outputParticlesInitial.empty() and
      not m_cfg.outputParticlesFinal.empty()) {
    // Initial state of particles
    // Register to the event store
    m_outputParticlesInitial(
        ctx, SimParticleContainer(eventData.particlesInitial.begin(),
                                  eventData.particlesInitial.end()));
    // Final state of particles
    // Register to the event store
    m_outputParticlesFinal(
        ctx, SimParticleContainer(eventData.particlesFinal.begin(),
                                  eventData.particlesFinal.end()));
  }

  // Output handling: Simulated hits
  if (not m_cfg.outputSimHits.empty()) {
    // Register to the event store
    m_outputSimHits(
        ctx, SimHitContainer(eventData.hits.begin(), eventData.hits.end()));
  }

  // Output handling: Material tracks
  if (not m_cfg.outputMaterialTracks.empty()) {
    m_outputMaterialTracks(
        ctx, decltype(eventData.materialTracks)(eventData.materialTracks));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
