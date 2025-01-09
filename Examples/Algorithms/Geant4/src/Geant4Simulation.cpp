// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/FpeMonitoring/FpeMonitor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleKillAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Geant4/SteppingActionList.hpp"

#include <stdexcept>
#include <utility>

#include <G4FieldManager.hh>
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
#include <Randomize.hh>
#include <boost/version.hpp>

namespace ActsExamples {

Geant4SimulationBase::Geant4SimulationBase(const Config& cfg, std::string name,
                                           Acts::Logging::Level level)
    : IAlgorithm(std::move(name), level) {
  if (cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (cfg.detector == nullptr) {
    throw std::invalid_argument("Missing detector construction factory");
  }
  if (cfg.randomNumbers == nullptr) {
    throw std::invalid_argument("Missing random numbers");
  }

  m_logger = Acts::getDefaultLogger("Geant4", level);

  m_eventStore = std::make_shared<Geant4::EventStore>();

  // tweak logging
  // If we are in VERBOSE mode, set the verbose level in Geant4 to 2.
  // 3 would be also possible, but that produces infinite amount of output.
  m_geant4Level = logger().level() == Acts::Logging::VERBOSE ? 2 : 0;
}

Geant4SimulationBase::~Geant4SimulationBase() = default;

void Geant4SimulationBase::commonInitialization() {
  // Set the detector construction
  {
    // Clear detector construction if it exists
    if (runManager().GetUserDetectorConstruction() != nullptr) {
      delete runManager().GetUserDetectorConstruction();
    }
    // G4RunManager will take care of deletion
    m_detectorConstruction =
        config()
            .detector
            ->buildGeant4DetectorConstruction(config().constructionOptions)
            .release();
    runManager().SetUserInitialization(m_detectorConstruction);
    runManager().InitializeGeometry();
  }

  m_geant4Instance->tweakLogging(m_geant4Level);
}

G4RunManager& Geant4SimulationBase::runManager() const {
  return *m_geant4Instance->runManager.get();
}

Geant4::EventStore& Geant4SimulationBase::eventStore() const {
  return *m_eventStore;
}

ProcessCode Geant4SimulationBase::initialize() {
  // Initialize the Geant4 run manager
  runManager().Initialize();

  return ProcessCode::SUCCESS;
}

ProcessCode Geant4SimulationBase::execute(const AlgorithmContext& ctx) const {
  // Ensure exclusive access to the Geant4 run manager
  std::lock_guard<std::mutex> guard(m_geant4Instance->mutex);

  // Set the seed new per event, so that we get reproducible results
  G4Random::setTheSeed(config().randomNumbers->generateSeed(ctx));

  // Get and reset event registry state
  eventStore() = Geant4::EventStore{};

  // Register the current event store to the registry
  // this will allow access from the User*Actions
  eventStore().store = &(ctx.eventStore);

  // Register the input particle read handle
  eventStore().inputParticles = &m_inputParticles;

  ACTS_DEBUG("Sending Geant RunManager the BeamOn() command.");
  {
    Acts::FpeMonitor mon{0};  // disable all FPEs while we're in Geant4
    // Start simulation. each track is simulated as a separate Geant4 event.
    runManager().BeamOn(1);
  }

  // Print out warnings about possible particle collision if happened
  if (eventStore().particleIdCollisionsInitial > 0 ||
      eventStore().particleIdCollisionsFinal > 0 ||
      eventStore().parentIdNotFound > 0) {
    ACTS_WARNING(
        "Particle ID collisions detected, don't trust the particle "
        "identification!");
    ACTS_WARNING(
        "- initial states: " << eventStore().particleIdCollisionsInitial);
    ACTS_WARNING("- final states: " << eventStore().particleIdCollisionsFinal);
    ACTS_WARNING("- parent ID not found: " << eventStore().parentIdNotFound);
  }

  if (eventStore().hits.empty()) {
    ACTS_DEBUG("Step merging: No steps recorded");
  } else {
    ACTS_DEBUG("Step merging: mean hits per hit: "
               << static_cast<double>(eventStore().numberGeantSteps) /
                      eventStore().hits.size());
    ACTS_DEBUG(
        "Step merging: max hits per hit: " << eventStore().maxStepsForHit);
  }

  return ProcessCode::SUCCESS;
}

std::shared_ptr<Geant4Handle> Geant4SimulationBase::geant4Handle() const {
  return m_geant4Instance;
}

Geant4Simulation::Geant4Simulation(const Config& cfg,
                                   Acts::Logging::Level level)
    : Geant4SimulationBase(cfg, "Geant4Simulation", level), m_cfg(cfg) {
  m_geant4Instance =
      m_cfg.geant4Handle
          ? m_cfg.geant4Handle
          : Geant4Manager::instance().createHandle(m_cfg.physicsList);
  if (m_geant4Instance->physicsListName != m_cfg.physicsList) {
    throw std::runtime_error("inconsistent physics list");
  }

  commonInitialization();

  // Set the primarty generator
  {
    // Clear primary generation action if it exists
    if (runManager().GetUserPrimaryGeneratorAction() != nullptr) {
      delete runManager().GetUserPrimaryGeneratorAction();
    }
    Geant4::SimParticleTranslation::Config prCfg;
    prCfg.eventStore = m_eventStore;
    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new Geant4::SimParticleTranslation(
        prCfg, m_logger->cloneWithSuffix("SimParticleTranslation"));
    // Set the primary generator action
    runManager().SetUserAction(primaryGeneratorAction);
  }

  // Particle action
  {
    // Clear tracking action if it exists
    if (runManager().GetUserTrackingAction() != nullptr) {
      delete runManager().GetUserTrackingAction();
    }
    Geant4::ParticleTrackingAction::Config trackingCfg;
    trackingCfg.eventStore = m_eventStore;
    trackingCfg.keepParticlesWithoutHits = cfg.keepParticlesWithoutHits;
    // G4RunManager will take care of deletion
    auto trackingAction = new Geant4::ParticleTrackingAction(
        trackingCfg, m_logger->cloneWithSuffix("ParticleTracking"));
    runManager().SetUserAction(trackingAction);
  }

  // Stepping actions
  Geant4::SensitiveSteppingAction* sensitiveSteppingActionAccess = nullptr;
  {
    // Clear stepping action if it exists
    if (runManager().GetUserSteppingAction() != nullptr) {
      delete runManager().GetUserSteppingAction();
    }

    Geant4::ParticleKillAction::Config particleKillCfg;
    particleKillCfg.eventStore = m_eventStore;
    particleKillCfg.volume = cfg.killVolume;
    particleKillCfg.maxTime = cfg.killAfterTime;
    particleKillCfg.secondaries = cfg.killSecondaries;

    Geant4::SensitiveSteppingAction::Config stepCfg;
    stepCfg.eventStore = m_eventStore;
    stepCfg.charged = cfg.recordHitsOfCharged;
    stepCfg.neutral = cfg.recordHitsOfNeutrals;
    stepCfg.primary = cfg.recordHitsOfPrimaries;
    stepCfg.secondary = cfg.recordHitsOfSecondaries;
    stepCfg.stepLogging = cfg.recordPropagationSummaries;

    Geant4::SteppingActionList::Config steppingCfg;
    steppingCfg.actions.push_back(std::make_unique<Geant4::ParticleKillAction>(
        particleKillCfg, m_logger->cloneWithSuffix("Killer")));

    auto sensitiveSteppingAction =
        std::make_unique<Geant4::SensitiveSteppingAction>(
            stepCfg, m_logger->cloneWithSuffix("SensitiveStepping"));
    sensitiveSteppingActionAccess = sensitiveSteppingAction.get();

    steppingCfg.actions.push_back(std::move(sensitiveSteppingAction));

    // G4RunManager will take care of deletion
    auto steppingAction = new Geant4::SteppingActionList(steppingCfg);
    runManager().SetUserAction(steppingAction);
  }

  // Get the g4World cache
  G4VPhysicalVolume* g4World = m_detectorConstruction->Construct();

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Set the magnetic field
  if (cfg.magneticField) {
    ACTS_INFO("Setting ACTS configured field to Geant4.");

    Geant4::MagneticFieldWrapper::Config g4FieldCfg;
    g4FieldCfg.magneticField = cfg.magneticField;
    m_magneticField =
        std::make_unique<Geant4::MagneticFieldWrapper>(g4FieldCfg);

    // Set the field or the G4Field manager
    m_fieldManager = std::make_unique<G4FieldManager>();
    m_fieldManager->SetDetectorField(m_magneticField.get());
    m_fieldManager->CreateChordFinder(m_magneticField.get());

    // Propagate down to all childrend
    g4World->GetLogicalVolume()->SetFieldManager(m_fieldManager.get(), true);
  }

  // ACTS sensitive surfaces are provided, so hit creation is turned on
  if (cfg.sensitiveSurfaceMapper != nullptr) {
    Geant4::SensitiveSurfaceMapper::State sState;
    ACTS_INFO(
        "Remapping selected volumes from Geant4 to Acts::Surface::GeometryID");
    cfg.sensitiveSurfaceMapper->remapSensitiveNames(
        sState, Acts::GeometryContext{}, g4World, Acts::Transform3::Identity());

    auto allSurfacesMapped = cfg.sensitiveSurfaceMapper->checkMapping(
        sState, Acts::GeometryContext{}, false, false);
    if (!allSurfacesMapped) {
      ACTS_WARNING(
          "Not all sensitive surfaces have been mapped to Geant4 volumes!");
    }

    sensitiveSteppingActionAccess->assignSurfaceMapping(
        sState.g4VolumeToSurfaces);
  }

  m_inputParticles.initialize(cfg.inputParticles);
  m_outputSimHits.initialize(cfg.outputSimHits);
  m_outputParticles.initialize(cfg.outputParticles);

  if (cfg.recordPropagationSummaries) {
    m_outputPropagationSummaries.initialize(cfg.outputPropagationSummaries);
  }
}

Geant4Simulation::~Geant4Simulation() = default;

ProcessCode Geant4Simulation::execute(const AlgorithmContext& ctx) const {
  auto ret = Geant4SimulationBase::execute(ctx);
  if (ret != ProcessCode::SUCCESS) {
    return ret;
  }

  // Output handling: Simulation
  m_outputParticles(
      ctx, SimParticleContainer(eventStore().particlesSimulated.begin(),
                                eventStore().particlesSimulated.end()));

#if BOOST_VERSION < 107800
  SimHitContainer container;
  for (const auto& hit : eventStore().hits) {
    container.insert(hit);
  }
  m_outputSimHits(ctx, std::move(container));
#else
  m_outputSimHits(
      ctx, SimHitContainer(eventStore().hits.begin(), eventStore().hits.end()));
#endif

  // Output the propagation summaries if requested
  if (m_cfg.recordPropagationSummaries) {
    PropagationSummaries summaries;
    summaries.reserve(eventStore().propagationRecords.size());
    for (auto& [trackId, summary] : eventStore().propagationRecords) {
      summaries.push_back(std::move(summary));
    }
    m_outputPropagationSummaries(ctx, std::move(summaries));
  }

  return ProcessCode::SUCCESS;
}

Geant4MaterialRecording::Geant4MaterialRecording(const Config& cfg,
                                                 Acts::Logging::Level level)
    : Geant4SimulationBase(cfg, "Geant4Simulation", level), m_cfg(cfg) {
  auto physicsListName = "MaterialPhysicsList";
  m_geant4Instance =
      m_cfg.geant4Handle
          ? m_cfg.geant4Handle
          : Geant4Manager::instance().createHandle(
                std::make_unique<Geant4::MaterialPhysicsList>(
                    m_logger->cloneWithSuffix("MaterialPhysicsList")),
                physicsListName);
  if (m_geant4Instance->physicsListName != physicsListName) {
    throw std::runtime_error("inconsistent physics list");
  }

  commonInitialization();

  // Set the primarty generator
  {
    // Clear primary generation action if it exists
    if (runManager().GetUserPrimaryGeneratorAction() != nullptr) {
      delete runManager().GetUserPrimaryGeneratorAction();
    }

    Geant4::SimParticleTranslation::Config prCfg;
    prCfg.eventStore = m_eventStore;
    prCfg.forcedPdgCode = 0;
    prCfg.forcedCharge = 0.;
    prCfg.forcedMass = 0.;

    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new Geant4::SimParticleTranslation(
        prCfg, m_logger->cloneWithSuffix("SimParticleTranslation"));
    // Set the primary generator action
    runManager().SetUserAction(primaryGeneratorAction);
  }

  // Particle action
  {
    // Clear tracking action if it exists
    if (runManager().GetUserTrackingAction() != nullptr) {
      delete runManager().GetUserTrackingAction();
    }
    Geant4::ParticleTrackingAction::Config trackingCfg;
    trackingCfg.eventStore = m_eventStore;
    trackingCfg.keepParticlesWithoutHits = true;
    // G4RunManager will take care of deletion
    auto trackingAction = new Geant4::ParticleTrackingAction(
        trackingCfg, m_logger->cloneWithSuffix("ParticleTracking"));
    runManager().SetUserAction(trackingAction);
  }

  // Stepping action
  {
    // Clear stepping action if it exists
    if (runManager().GetUserSteppingAction() != nullptr) {
      delete runManager().GetUserSteppingAction();
    }
    Geant4::MaterialSteppingAction::Config steppingCfg;
    steppingCfg.eventStore = m_eventStore;
    steppingCfg.excludeMaterials = m_cfg.excludeMaterials;
    // G4RunManager will take care of deletion
    auto steppingAction = new Geant4::MaterialSteppingAction(
        steppingCfg, m_logger->cloneWithSuffix("MaterialSteppingAction"));
    runManager().SetUserAction(steppingAction);
  }

  runManager().Initialize();

  m_inputParticles.initialize(cfg.inputParticles);
  m_outputMaterialTracks.initialize(cfg.outputMaterialTracks);
}

Geant4MaterialRecording::~Geant4MaterialRecording() = default;

ProcessCode Geant4MaterialRecording::execute(
    const AlgorithmContext& ctx) const {
  const auto ret = Geant4SimulationBase::execute(ctx);
  if (ret != ProcessCode::SUCCESS) {
    return ret;
  }

  // Output handling: Material tracks
  m_outputMaterialTracks(
      ctx, decltype(eventStore().materialTracks)(eventStore().materialTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
