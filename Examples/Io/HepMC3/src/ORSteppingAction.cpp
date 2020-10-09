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

#include <HepMC3/Attribute.h>

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
	//~ ParticleRecord p;
	//~ p.position[0] = step->GetPostStepPoint()->GetPosition().x() / CLHEP::mm;
	//~ p.position[1] = step->GetPostStepPoint()->GetPosition().y() / CLHEP::mm;
	//~ p.position[2] = step->GetPostStepPoint()->GetPosition().z() / CLHEP::mm;
	//~ p.momentum[0] = step->GetPostStepPoint()->GetMomentum().x() / CLHEP::GeV;
	//~ p.momentum[1] = step->GetPostStepPoint()->GetMomentum().y() / CLHEP::GeV;
	//~ p.momentum[2] = step->GetPostStepPoint()->GetMomentum().z() / CLHEP::GeV;
	//~ p.globalTime = step->GetPostStepPoint()->GetGlobalTime();
	//~ p.pdg = step->GetTrack()->GetDynamicParticle()->GetPDGcode();
	//~ p.vertex[0] = step->GetTrack()->GetVertexPosition().x();
	//~ p.vertex[1] = step->GetTrack()->GetVertexPosition().y();
	//~ p.vertex[2] = step->GetTrack()->GetVertexPosition().z();
	//~ p.energy = step->GetPostStepPoint()->GetTotalEnergy() / CLHEP::GeV;
	//~ p.mass = step->GetPostStepPoint()->GetMass();
	//~ p.charge = step->GetPostStepPoint()->GetCharge();
	//~ p.trackid = step->GetTrack()->GetTrackID();
	//~ p.parentid = step->GetTrack()->GetParentID();
	//~ p.volume = (step->GetPostStepPoint()->GetPhysicalVolume() != nullptr) ? step->GetPostStepPoint()->GetPhysicalVolume()->GetName() : "No volume";
	//~ p.process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	
	auto* track = step->GetTrack();

	// The particle after the step
	auto postStepMomentum = track->GetMomentum();
	auto postStepEnergy = track->GetTotalEnergy();
	HepMC3::FourVector mom4{postStepMomentum[0], postStepMomentum[1], postStepMomentum[2], postStepEnergy};
	auto postParticle = std::make_shared<HepMC3::GenParticle>(mom4, track->GetDynamicParticle()->GetPDGcode());
	
	auto process = std::make_shared<HepMC3::StringAttribute>(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
	
	// Assign particle to production vertex
	if(!m_previousVertex)
	{
		// Get the production position of the particle
		auto* preStep = step->GetPreStepPoint();
		auto prePosition = preStep->GetPosition();
		auto preTime = preStep->GetGlobalTime();
		HepMC3::FourVector prePos{prePosition[0], prePosition[1], prePosition[2], preTime};
	
		// Handle the first step
		if(OREventAction::instance()->event()->vertices().empty())
		{
			auto vertex = std::make_shared<HepMC3::GenVertex>(prePos);
			vertex->add_particle_out(postParticle);
			OREventAction::instance()->event()->add_vertex(vertex);
			vertex->add_attribute("Process-" + std::to_string(postParticle->id()), process);
		}
		else
			// Search for an existing vertex
			for(const auto& vertex : OREventAction::instance()->event()->vertices())
			{
				if(vertex->position() == prePos)
				{
					vertex->add_particle_out(postParticle); // TODO: the initial parameters need to be stored additionally
					vertex->add_attribute("Process-" + std::to_string(postParticle->id()), process);
				}
			}
	}
	else
	{
		// Add particle from same track to vertex
		m_previousVertex->add_particle_out(postParticle);
		m_previousVertex->add_attribute("Process-" + std::to_string(postParticle->id()), process);
	}

	// Build the vertex after this step
	auto* postStep = step->GetPostStepPoint();	
	auto postPosition = postStep->GetPosition();
	auto postTime = postStep->GetGlobalTime();
	HepMC3::FourVector postPos{postPosition[0], postPosition[1], postPosition[2], postTime};
	m_previousVertex = std::make_shared<HepMC3::GenVertex>(postPos);
	
	// Add properties to the vertex
	m_previousVertex->add_particle_in(postParticle);
	
	// Store the vertex
	OREventAction::instance()->event()->add_vertex(m_previousVertex);
	
	// Stop tracking the vertex if the particle dies
	if(track->GetTrackStatus() != fAlive)
	{
		process = std::make_shared<HepMC3::StringAttribute>("Death");
		m_previousVertex->add_attribute("Process-" + std::to_string(postParticle->id()), process);
		m_previousVertex = nullptr;
	}
	
testCounter++;	
//~ std::cout << "energy: " << track->GetKineticEnergy() << " | " << track->GetVertexKineticEnergy() << std::endl;
//~ std::cout << track->GetParentID() << " | " << track->GetTrackID() << " | " << track->GetDynamicParticle()->GetPDGcode() 
    //~ << " | " << track->GetTrackLength() << " | " << step->GetPreStepPoint()->GetLocalTime() 
	//~ << " | " << (step->GetPreStepPoint()->GetProcessDefinedStep() ? step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() : "---") 
	//~ << " | " << (track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "---") << " | " << track->GetCurrentStepNumber() << std::endl;	
}

void
ActsExamples::ORSteppingAction::clear()
{
	m_previousVertex = nullptr;
}