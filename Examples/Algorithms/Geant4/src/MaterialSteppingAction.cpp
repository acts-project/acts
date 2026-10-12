// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"

#include "Acts/Material/MaterialInteraction.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"

#include <cassert>
#include <cstddef>
#include <optional>
#include <ostream>
#include <unordered_map>
#include <utility>

#include <G4Step.hh>
#include <G4Track.hh>

namespace ActsExamples::Geant4 {

MaterialSteppingAction::MaterialSteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(),
      m_cfg(cfg),
      m_stepConverter(cfg.stepConverter),
      m_logger(std::move(logger)) {}

void MaterialSteppingAction::UserSteppingAction(const G4Step* stepPtr) {
  assert(stepPtr != nullptr);
  const G4Step& step = *stepPtr;

  std::optional<Acts::MaterialInteraction> mInteraction =
      m_stepConverter.convert(step);
  if (!mInteraction.has_value()) {
    ACTS_VERBOSE("Skip step without recorded material.");
    return;
  }

  ACTS_VERBOSE(
      "Performing a step with step size = " << mInteraction->pathCorrection);

  assert(step.GetTrack() != nullptr);
  const G4Track& track = *step.GetTrack();
  const std::size_t trackId = track.GetTrackID();
  auto [it, inserted] = eventStore().materialTracks.try_emplace(trackId - 1);
  if (inserted) {
    it->second.first =
        ActsPlugins::Geant4MaterialStepConverter::trackStart(track);
  }
  it->second.second.materialInteractions.push_back(std::move(*mInteraction));
}

}  // namespace ActsExamples::Geant4
