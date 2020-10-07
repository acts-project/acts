// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ORSteppingAction.hpp"
#include "OREventAction.hpp"
#include <stdexcept>
#include "Acts/Utilities/Units.hpp"
#include "G4Step.hh"
#include "G4VProcess.hh"
//~ #include "SystemOfUnits.h"

ActsExamples::ORSteppingAction* ActsExamples::ORSteppingAction::s_instance
    = nullptr;

ActsExamples::ORSteppingAction*
ActsExamples::ORSteppingAction::instance()
{
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::ORSteppingAction::ORSteppingAction()
  : G4UserSteppingAction()
{
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::ORSteppingAction::~ORSteppingAction()
{
  s_instance = nullptr;
}

void
ActsExamples::ORSteppingAction::UserSteppingAction(const G4Step* step)
{
	ParticleRecord p;
	p.position[0] = step->GetPostStepPoint()->GetPosition().x() / CLHEP::mm;
	p.position[1] = step->GetPostStepPoint()->GetPosition().y() / CLHEP::mm;
	p.position[2] = step->GetPostStepPoint()->GetPosition().z() / CLHEP::mm;
	p.momentum[0] = step->GetPostStepPoint()->GetMomentum().x() / CLHEP::GeV;
	p.momentum[1] = step->GetPostStepPoint()->GetMomentum().y() / CLHEP::GeV;
	p.momentum[2] = step->GetPostStepPoint()->GetMomentum().z() / CLHEP::GeV;
	p.globalTime = step->GetPostStepPoint()->GetGlobalTime();
	p.pdg = step->GetTrack()->GetDynamicParticle()->GetPDGcode();
	p.vertex[0] = step->GetTrack()->GetVertexPosition().x();
	p.vertex[1] = step->GetTrack()->GetVertexPosition().y();
	p.vertex[2] = step->GetTrack()->GetVertexPosition().z();
	p.energy = step->GetPostStepPoint()->GetTotalEnergy() / CLHEP::GeV;
	p.mass = step->GetPostStepPoint()->GetMass();
	p.charge = step->GetPostStepPoint()->GetCharge();
	p.trackid = step->GetTrack()->GetTrackID();
	p.parentid = step->GetTrack()->GetParentID();
	p.volume = (step->GetPostStepPoint()->GetPhysicalVolume() != nullptr) ? step->GetPostStepPoint()->GetPhysicalVolume()->GetName() : "No volume";
	p.process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
//~ std::cout << p.process << std::endl;	
	
	m_particles[p.trackid].push_back(p);
}

void
ActsExamples::ORSteppingAction::clear()
{
	m_particles.clear();
}