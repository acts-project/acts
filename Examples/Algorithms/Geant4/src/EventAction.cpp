// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"

#include "ActsExamples/Geant4/PrimaryGeneratorAction.hpp"

#include <stdexcept>

#include <G4Event.hh>
#include <G4RunManager.hh>

#include "SteppingAction.hpp"

namespace ActsExamples::Geant4 {

EventAction* EventAction::s_instance = nullptr;

EventAction* EventAction::instance() {
  return s_instance;
}

EventAction::EventAction() : G4UserEventAction() {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate the EventAction singleton");
  } else {
    s_instance = this;
  }
}

EventAction::~EventAction() {
  s_instance = nullptr;
}

void EventAction::BeginOfEventAction(const G4Event*) {
  // reset the collection of material steps
  SteppingAction::instance()->clear();
}

void EventAction::EndOfEventAction(const G4Event* event) {
  const auto* rawPos = event->GetPrimaryVertex();
  // access the initial direction of the track
  G4ThreeVector rawDir = PrimaryGeneratorAction::instance()->direction();
  // create the RecordedMaterialTrack
  Acts::RecordedMaterialTrack mtrecord;
  mtrecord.first.first =
      Acts::Vector3(rawPos->GetX0(), rawPos->GetY0(), rawPos->GetZ0());
  mtrecord.first.second = Acts::Vector3(rawDir.x(), rawDir.y(), rawDir.z());
  mtrecord.second.materialInteractions =
      SteppingAction::instance()->materialSteps();

  // write out the RecordedMaterialTrack of one event
  m_materialTracks.push_back(mtrecord);
}

/// Clear the recorded data.
void EventAction::clear() {
  m_materialTracks.clear();
}

/// Access the recorded material tracks.
///
/// This only contains valid data after the end-of-event action has been
/// executed.
const std::vector<Acts::RecordedMaterialTrack>& EventAction::materialTracks()
    const {
  return m_materialTracks;
}

}  // namespace ActsExamples::Geant4