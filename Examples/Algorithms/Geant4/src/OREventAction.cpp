// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/OREventAction.hpp"
#include <stdexcept>
#include "ActsExamples/Geant4/ORPrimaryGeneratorAction.hpp"
#include "ActsExamples/Geant4/ORSteppingAction.hpp"
#include <G4Event.hh>
#include <G4RunManager.hh>

ActsExamples::OREventAction* ActsExamples::OREventAction::s_instance = nullptr;

ActsExamples::OREventAction*
ActsExamples::OREventAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

ActsExamples::OREventAction::OREventAction() : G4UserEventAction()
{
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::OREventAction::~OREventAction()
{
  s_instance = nullptr;
}

void
ActsExamples::OREventAction::BeginOfEventAction(const G4Event*)
{
  // reset the collection of material steps
  //~ m_particles.clear(); // TODO: Is this one required?
  ORSteppingAction::Instance()->Reset();
}

void
ActsExamples::OREventAction::EndOfEventAction(const G4Event*)
{
	// TODO: Wrap up tracks to build collection
}

void
ActsExamples::OREventAction::clear()
{
	m_particles.clear();
}

//~ void EventAction::EndOfEventAction(const G4Event* event) {
  //~ const auto* rawPos = event->GetPrimaryVertex();
  //~ // access the initial direction of the track
  //~ G4ThreeVector rawDir = PrimaryGeneratorAction::instance()->direction();
  //~ // create the RecordedMaterialTrack
  //~ Acts::RecordedMaterialTrack mtrecord;
  //~ mtrecord.first.first =
      //~ Acts::Vector3D(rawPos->GetX0(), rawPos->GetY0(), rawPos->GetZ0());
  //~ mtrecord.first.second = Acts::Vector3D(rawDir.x(), rawDir.y(), rawDir.z());
  //~ mtrecord.second.materialInteractions =
      //~ SteppingAction::instance()->materialSteps();

  //~ // write out the RecordedMaterialTrack of one event
  //~ m_materialTracks.push_back(mtrecord);
//~ }

/// Access the recorded material tracks.
///
/// This only contains valid data after the end-of-event action has been
/// executed.
const std::vector<Acts::RecordedMaterialTrack>& EventAction::materialTracks()
    const {
  return m_materialTracks;
}
