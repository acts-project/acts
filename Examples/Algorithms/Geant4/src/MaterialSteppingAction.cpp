// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"

#include "ActsExamples/Geant4/EventStore.hpp"

#include <cassert>
#include <cstddef>
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
      m_stepConverter({.excludeMaterials = cfg.excludeMaterials,
                       .recordElementFractions = cfg.recordElementFractions}),
      m_logger(std::move(logger)) {}

void MaterialSteppingAction::UserSteppingAction(const G4Step* stepPtr) {
  assert(stepPtr != nullptr);
  const G4Step& step = *stepPtr;
  assert(step.GetTrack() != nullptr);
  const G4Track& track = *step.GetTrack();

  // Every track gets an entry, also if it never crosses recorded material
  auto& materialTracks = eventStore().materialTracks;
  const std::size_t trackId = track.GetTrackID() - 1;
  auto it = materialTracks.find(trackId);
  if (it == materialTracks.end()) {
    it = materialTracks.emplace(trackId, m_stepConverter.startTrack(track))
             .first;
  }

  if (m_stepConverter.addStep(step, it->second)) {
    ACTS_VERBOSE("Recorded a step with step size = " << step.GetStepLength());
  } else {
    ACTS_VERBOSE("Skip step without recorded material.");
  }
}

}  // namespace ActsExamples::Geant4
