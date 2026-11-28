// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "SimulationActions.hpp"

#include "ActsExamples/Geant4/ParticleKillAction.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples::Geant4 {

SimulationActions::SimulationActions(
    const Geant4Simulation::Config& cfg,
    const std::shared_ptr<EventContext>& context, const Acts::Logger& logger) {
  SimParticleTranslation::Config primaryCfg;
  primaryCfg.eventContext = context;
  primaryCfg.preserveGeant4EventId = cfg.numGeant4Workers > 0;
  generator = std::make_unique<SimParticleTranslation>(
      primaryCfg, logger.cloneWithSuffix("SimParticleTranslation"));

  ParticleTrackingAction::Config trackingCfg;
  trackingCfg.eventContext = context;
  trackingCfg.keepParticlesWithoutHits = cfg.keepParticlesWithoutHits;
  tracking = std::make_unique<ParticleTrackingAction>(
      trackingCfg, logger.cloneWithSuffix("ParticleTracking"));

  ParticleKillAction::Config killCfg;
  killCfg.eventContext = context;
  killCfg.volume = cfg.killVolume;
  killCfg.maxTime = cfg.killAfterTime;
  killCfg.secondaries = cfg.killSecondaries;

  SensitiveSteppingAction::Config hitCfg;
  hitCfg.eventContext = context;
  hitCfg.charged = cfg.recordHitsOfCharged;
  hitCfg.neutral = cfg.recordHitsOfNeutrals;
  hitCfg.primary = cfg.recordHitsOfPrimaries;
  hitCfg.secondary = cfg.recordHitsOfSecondaries;
  hitCfg.stepLogging = cfg.recordPropagationSummaries;
  auto hits = std::make_shared<SensitiveSteppingAction>(
      hitCfg, logger.cloneWithSuffix("SensitiveStepping"));
  sensitive = hits.get();
  SteppingActionList::Config steps;
  steps.actions.push_back(std::make_shared<ParticleKillAction>(
      killCfg, logger.cloneWithSuffix("Killer")));
  steps.actions.push_back(std::move(hits));
  stepping = std::make_unique<SteppingActionList>(steps);
}

void appendDiagnostics(SteppingActionList& actions,
                       const Geant4Simulation::Config& cfg,
                       const std::shared_ptr<EventContext>& context) {
  SteppingActionList::Config extensions;
  for (const auto& factory : cfg.steppingActionFactories) {
    auto action = factory(context);
    if (!action) {
      throw std::invalid_argument(
          "Geant4 stepping action factory returned null");
    }
    extensions.actions.push_back(std::move(action));
  }
  actions.appendActions(extensions);
}

SensitiveSurfaceMapper::VolumeToSurfAssocMap_t mapSensitiveVolumes(
    const Geant4Simulation::Config& cfg, G4VPhysicalVolume& world,
    const Acts::Logger& logger) {
  SensitiveSurfaceMapper::State state;
  if (cfg.sensitiveSurfaceMapper) {
    auto geometry = Acts::GeometryContext::dangerouslyDefaultConstruct();
    cfg.sensitiveSurfaceMapper->remapSensitiveNames(
        state, geometry, &world, Acts::Transform3::Identity());
    if (!cfg.sensitiveSurfaceMapper->checkMapping(state, geometry, false,
                                                  false)) {
      ACTS_LOG_WITH_LOGGER(
          logger, Acts::Logging::WARNING,
          "Not all sensitive surfaces have been mapped to Geant4 volumes!");
    }
  }
  return std::move(state.g4VolumeToSurfaces);
}

}  // namespace ActsExamples::Geant4
