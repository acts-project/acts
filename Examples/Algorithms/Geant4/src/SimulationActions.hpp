// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Geant4/SteppingActionList.hpp"

#include <memory>

namespace ActsExamples::Geant4 {

/// Actions for one serial simulation instance or one native transport worker.
struct SimulationActions {
  SimulationActions(const Geant4Simulation::Config& cfg,
                    const std::shared_ptr<EventContext>& context,
                    const Acts::Logger& logger);

  std::unique_ptr<SimParticleTranslation> generator;
  std::unique_ptr<ParticleTrackingAction> tracking;
  std::unique_ptr<SteppingActionList> stepping;
  // Owned by the stepping list.
  SensitiveSteppingAction* sensitive;
};

void appendDiagnostics(SteppingActionList& actions,
                       const Geant4Simulation::Config& cfg,
                       const std::shared_ptr<EventContext>& context);

SensitiveSurfaceMapper::VolumeToSurfAssocMap_t mapSensitiveVolumes(
    const Geant4Simulation::Config& cfg, G4VPhysicalVolume& world,
    const Acts::Logger& logger);

}  // namespace ActsExamples::Geant4
