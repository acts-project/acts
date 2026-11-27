// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
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

#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserLimits.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4Version.hh>

#include "SimulationSession.hpp"

namespace ActsExamples {

Geant4SimulationBase::Geant4SimulationBase(
    const Config& cfg, const std::string& name,
    std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm(name, std::move(logger)) {
  if (cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (cfg.detector == nullptr) {
    throw std::invalid_argument("Missing detector construction factory");
  }
  if (cfg.randomNumbers == nullptr) {
    throw std::invalid_argument("Missing random numbers");
  }

  // tweak logging
  // If we are in VERBOSE mode, set the verbose level in Geant4 to 2.
  // 3 would be also possible, but that produces infinite amount of output.
  m_geant4Level = this->logger().level() == Acts::Logging::VERBOSE ? 2 : 0;
}

Geant4SimulationBase::~Geant4SimulationBase() {
  // Sequencer finalization is skipped if an earlier algorithm throws.
  try {
    finalize();
  } catch (const std::exception& e) {
    try {
      ACTS_ERROR("Failed to close Geant4 run during destruction: " << e.what());
    } catch (...) {
      // Logging can throw; destruction must preserve the original failure.
    }
  } catch (...) {
    // A user action can throw a non-standard exception during cleanup.
  }
}

void Geant4SimulationBase::commonInitialization() {
  auto construction = config().detector->buildGeant4DetectorConstruction(
      config().constructionOptions);
  m_session->install();
  const auto* previous = runManager().GetUserDetectorConstruction();
  m_detectorConstruction = construction.release();
  runManager().SetUserInitialization(m_detectorConstruction);
  delete previous;
  runManager().InitializeGeometry();

  geant4Handle()->tweakLogging(m_geant4Level);
}

G4RunManager& Geant4SimulationBase::runManager() const {
  return *geant4Handle()->runManager;
}

std::shared_ptr<Geant4::EventContext> Geant4SimulationBase::eventContext()
    const {
  return m_session->eventContext;
}

ProcessCode Geant4SimulationBase::initialize() {
  m_session->beginRun();
  return ProcessCode::SUCCESS;
}

ProcessCode Geant4SimulationBase::finalize() {
  if (m_session) {
    m_session->endRun();
  }
  return ProcessCode::SUCCESS;
}

Geant4::EventResult Geant4SimulationBase::simulateEvent(
    const AlgorithmContext& ctx) const {
  ACTS_DEBUG("Processing Geant4 event");
  auto result = m_session->simulate(m_inputParticles(ctx), ctx.simGeoContext,
                                    config().randomNumbers->generateSeed(ctx));

  // Print out warnings about possible particle collision if happened
  if (result.particleIdCollisionsInitial > 0 ||
      result.particleIdCollisionsFinal > 0 || result.parentIdNotFound > 0) {
    ACTS_WARNING(
        "Particle ID collisions detected, don't trust the particle "
        "identification!");
    ACTS_WARNING("- initial states: " << result.particleIdCollisionsInitial);
    ACTS_WARNING("- final states: " << result.particleIdCollisionsFinal);
    ACTS_WARNING("- parent ID not found: " << result.parentIdNotFound);
  }

  if (result.hits.empty()) {
    ACTS_DEBUG("Step merging: No steps recorded");
  } else {
    ACTS_DEBUG("Step merging: mean hits per hit: "
               << static_cast<double>(result.numberGeantSteps) /
                      result.hits.size());
    ACTS_DEBUG("Step merging: max hits per hit: " << result.maxStepsForHit);
  }

  return result;
}

std::shared_ptr<Geant4Handle> Geant4SimulationBase::geant4Handle() const {
  return m_session->handle;
}

Geant4Simulation::Geant4Simulation(const Config& cfg,
                                   std::unique_ptr<const Acts::Logger> logger)
    : Geant4SimulationBase(cfg, "Geant4Simulation", std::move(logger)),
      m_cfg(cfg) {
  for (const auto& factory : cfg.steppingActionFactories) {
    if (!factory) {
      throw std::invalid_argument("Empty Geant4 stepping action factory");
    }
  }
  m_inputParticles.initialize(cfg.inputParticles);
  m_outputSimHits.initialize(cfg.outputSimHits);
  m_outputParticles.initialize(cfg.outputParticles);

  if (cfg.recordPropagationSummaries) {
    m_outputPropagationSummaries.initialize(cfg.outputPropagationSummaries);
  }

  auto handle = m_cfg.geant4Handle
                    ? m_cfg.geant4Handle
                    : Geant4Manager::instance().createHandle(m_cfg.physicsList);
  m_session = std::make_unique<Geant4::SimulationSession>(
      std::move(handle), m_cfg.physicsList,
      cfg.propagatorLargestAcceptableStep);

  commonInitialization();

  // Set the primarty generator
  {
    Geant4::SimParticleTranslation::Config prCfg;
    prCfg.eventContext = eventContext();
    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new Geant4::SimParticleTranslation(
        prCfg, this->logger().cloneWithSuffix("SimParticleTranslation"));
    // Set the primary generator action
    const auto* previous = runManager().GetUserPrimaryGeneratorAction();
    runManager().SetUserAction(primaryGeneratorAction);
    delete previous;
  }

  // Particle action
  {
    Geant4::ParticleTrackingAction::Config trackingCfg;
    trackingCfg.eventContext = eventContext();
    trackingCfg.keepParticlesWithoutHits = cfg.keepParticlesWithoutHits;
    // G4RunManager will take care of deletion
    auto trackingAction = new Geant4::ParticleTrackingAction(
        trackingCfg, this->logger().cloneWithSuffix("ParticleTracking"));
    const auto* previous = runManager().GetUserTrackingAction();
    runManager().SetUserAction(trackingAction);
    delete previous;
  }

  // Stepping actions
  Geant4::SensitiveSteppingAction* sensitiveSteppingActionAccess = nullptr;
  {
    Geant4::ParticleKillAction::Config particleKillCfg;
    particleKillCfg.eventContext = eventContext();
    particleKillCfg.volume = cfg.killVolume;
    particleKillCfg.maxTime = cfg.killAfterTime;
    particleKillCfg.secondaries = cfg.killSecondaries;

    Geant4::SensitiveSteppingAction::Config stepCfg;
    stepCfg.eventContext = eventContext();
    stepCfg.charged = cfg.recordHitsOfCharged;
    stepCfg.neutral = cfg.recordHitsOfNeutrals;
    stepCfg.primary = cfg.recordHitsOfPrimaries;
    stepCfg.secondary = cfg.recordHitsOfSecondaries;
    stepCfg.stepLogging = cfg.recordPropagationSummaries;

    Geant4::SteppingActionList::Config steppingCfg;
    steppingCfg.actions.push_back(std::make_unique<Geant4::ParticleKillAction>(
        particleKillCfg, this->logger().cloneWithSuffix("Killer")));

    auto sensitiveSteppingAction =
        std::make_unique<Geant4::SensitiveSteppingAction>(
            stepCfg, this->logger().cloneWithSuffix("SensitiveStepping"));
    sensitiveSteppingActionAccess = sensitiveSteppingAction.get();

    steppingCfg.actions.push_back(std::move(sensitiveSteppingAction));
    // G4RunManager will take care of deletion
    auto steppingAction = new Geant4::SteppingActionList(steppingCfg);
    const auto* previous = runManager().GetUserSteppingAction();
    runManager().SetUserAction(steppingAction);
    delete previous;

    if (!cfg.steppingActionFactories.empty()) {
      m_session->configureActions = [this, steppingAction] {
        Geant4::SteppingActionList::Config extensions;
        for (const auto& factory : m_cfg.steppingActionFactories) {
          auto action = factory(eventContext());
          if (!action) {
            throw std::invalid_argument(
                "Geant4 stepping action factory returned null");
          }
          extensions.actions.push_back(std::move(action));
        }
        steppingAction->appendActions(extensions);
      };
    }
  }

  // Get the g4World cache
  G4VPhysicalVolume* g4World = m_detectorConstruction->Construct();

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Set the magnetic field
  if (cfg.magneticField) {
    ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::INFO,
                         "Setting ACTS configured field to Geant4.");

    m_session->setMagneticField(*g4World, cfg.magneticField);
  }

  // ACTS sensitive surfaces are provided, so hit creation is turned on
  if (cfg.sensitiveSurfaceMapper != nullptr) {
    Geant4::SensitiveSurfaceMapper::State sState;
    ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::INFO,
                         "Remapping selected volumes from Geant4 to "
                         "Acts::Surface::GeometryID");
    cfg.sensitiveSurfaceMapper->remapSensitiveNames(
        sState, Acts::GeometryContext::dangerouslyDefaultConstruct(), g4World,
        Acts::Transform3::Identity());

    auto allSurfacesMapped = cfg.sensitiveSurfaceMapper->checkMapping(
        sState, Acts::GeometryContext::dangerouslyDefaultConstruct(), false,
        false);
    if (!allSurfacesMapped) {
      ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::WARNING,
                           "Not all sensitive surfaces have been mapped to "
                           "Geant4 volumes!");
    }

    sensitiveSteppingActionAccess->assignSurfaceMapping(
        sState.g4VolumeToSurfaces);
  }
}

Geant4Simulation::~Geant4Simulation() = default;

ProcessCode Geant4Simulation::execute(const AlgorithmContext& ctx) const {
  auto result = simulateEvent(ctx);

  // Output handling: Simulation
  m_outputParticles(ctx, SimParticleContainer(result.particlesSimulated.begin(),
                                              result.particlesSimulated.end()));

  SimHitContainer hits;
  hits.adopt_sequence(std::move(result.hits));
  m_outputSimHits(ctx, std::move(hits));

  // Output the propagation summaries if requested
  if (m_cfg.recordPropagationSummaries) {
    PropagationSummaries summaries;
    summaries.reserve(result.propagationRecords.size());
    for (auto& [trackId, summary] : result.propagationRecords) {
      summaries.push_back(std::move(summary));
    }
    m_outputPropagationSummaries(ctx, std::move(summaries));
  }

  return ProcessCode::SUCCESS;
}

Geant4MaterialRecording::Geant4MaterialRecording(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : Geant4SimulationBase(cfg, "Geant4Simulation", std::move(logger)),
      m_cfg(cfg) {
  m_inputParticles.initialize(cfg.inputParticles);
  m_outputMaterialTracks.initialize(cfg.outputMaterialTracks);

  auto physicsListName = "MaterialPhysicsList";
  auto handle =
      m_cfg.geant4Handle
          ? m_cfg.geant4Handle
          : Geant4Manager::instance().createHandle(
                std::make_unique<Geant4::MaterialPhysicsList>(
                    this->logger().cloneWithSuffix("MaterialPhysicsList")),
                physicsListName);
  m_session = std::make_unique<Geant4::SimulationSession>(
      std::move(handle), physicsListName, cfg.propagatorLargestAcceptableStep);

  commonInitialization();

  // Set the primarty generator
  {
    Geant4::SimParticleTranslation::Config prCfg;
    prCfg.eventContext = eventContext();
    prCfg.forcedPdgCode = 0;
    prCfg.forcedCharge = 0.;
    prCfg.forcedMass = 0.;

    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new Geant4::SimParticleTranslation(
        prCfg, this->logger().cloneWithSuffix("SimParticleTranslation"));
    // Set the primary generator action
    const auto* previous = runManager().GetUserPrimaryGeneratorAction();
    runManager().SetUserAction(primaryGeneratorAction);
    delete previous;
  }

  // Particle action
  {
    Geant4::ParticleTrackingAction::Config trackingCfg;
    trackingCfg.eventContext = eventContext();
    trackingCfg.keepParticlesWithoutHits = true;
    // G4RunManager will take care of deletion
    auto trackingAction = new Geant4::ParticleTrackingAction(
        trackingCfg, this->logger().cloneWithSuffix("ParticleTracking"));
    const auto* previous = runManager().GetUserTrackingAction();
    runManager().SetUserAction(trackingAction);
    delete previous;
  }

  // Stepping action
  {
    Geant4::MaterialSteppingAction::Config steppingCfg;
    steppingCfg.eventContext = eventContext();
    steppingCfg.excludeMaterials = m_cfg.excludeMaterials;
    steppingCfg.recordElementFractions = m_cfg.recordElementFractions;
    // G4RunManager will take care of deletion
    auto steppingAction = new Geant4::MaterialSteppingAction(
        steppingCfg, this->logger().cloneWithSuffix("MaterialSteppingAction"));
    const auto* previous = runManager().GetUserSteppingAction();
    runManager().SetUserAction(steppingAction);
    delete previous;
  }

  runManager().Initialize();
}

Geant4MaterialRecording::~Geant4MaterialRecording() = default;

ProcessCode Geant4MaterialRecording::execute(
    const AlgorithmContext& ctx) const {
  auto result = simulateEvent(ctx);
  m_outputMaterialTracks(ctx, std::move(result.materialTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
