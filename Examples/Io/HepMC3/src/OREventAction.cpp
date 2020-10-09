// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "OREventAction.hpp"
#include <stdexcept>
#include "ORPrimaryGeneratorAction.hpp"
#include "ORSteppingAction.hpp"
#include <G4Event.hh>
#include <G4RunManager.hh>

ActsExamples::OREventAction* ActsExamples::OREventAction::s_instance = nullptr;

ActsExamples::OREventAction*
ActsExamples::OREventAction::instance()
{
  // Static acces function via G4RunManager
  return s_instance;
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
  ORSteppingAction::instance()->clear(); // TODO: will this remain?
  
  m_event = std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
}

void
ActsExamples::OREventAction::EndOfEventAction(const G4Event*)
{
	// TODO: Wrap up tracks to build collection
	std::cout << "Particles: " << m_event->particles().size() << " | " << "Vertices: " << m_event->vertices().size() << std::endl;
	std::cout << "Number of steps: " << ORSteppingAction::instance()->counter() << std::endl;
	ORSteppingAction::instance()->counter() = 0;
}

void
ActsExamples::OREventAction::clear()
{
	m_event = nullptr;
	// TODO | m_processTracks.clear();
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
