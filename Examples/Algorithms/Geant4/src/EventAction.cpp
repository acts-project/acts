// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"

#include <G4Event.hh>
#include <G4RunManager.hh>
#include <stdexcept>

#include "PrimaryGeneratorAction.hpp"
#include "SteppingAction.hpp"

using namespace ActsExamples;

EventAction* EventAction::s_instance = nullptr;

EventAction* EventAction::Instance() {
  // Static acces function via G4RunManager
  return s_instance;
}

EventAction::EventAction() : G4UserEventAction() {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

EventAction::~EventAction() {
  s_instance = nullptr;
}

void EventAction::BeginOfEventAction(const G4Event*) {
  // reset the collection of material steps
  SteppingAction::Instance()->Reset();
}

void EventAction::EndOfEventAction(const G4Event* event) {
  const auto* rawPos = event->GetPrimaryVertex();
  // access the initial direction of the track
  G4ThreeVector rawDir = PrimaryGeneratorAction::Instance()->direction();
  // create the RecordedMaterialTrack
  Acts::RecordedMaterialTrack mtrecord;
  mtrecord.first.first =
      Acts::Vector3D(rawPos->GetX0(), rawPos->GetY0(), rawPos->GetZ0());
  mtrecord.first.second = Acts::Vector3D(rawDir.x(), rawDir.y(), rawDir.z());
  mtrecord.second.materialInteractions =
      SteppingAction::Instance()->materialSteps();

  // write out the RecordedMaterialTrack of one event
  m_records.push_back(mtrecord);

  // write out the steps of one track in an event
  m_tracksteps.adopt_sequence(SteppingAction::Instance()->trackSteps());
}

void EventAction::Reset() {}
