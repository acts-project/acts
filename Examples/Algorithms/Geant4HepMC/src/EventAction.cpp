// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"
#include <stdexcept>
#include <G4Event.hh>
#include <G4RunManager.hh>
#include "SteppingAction.hpp"

ActsExamples::EventAction* ActsExamples::EventAction::s_instance = nullptr;

ActsExamples::EventAction* ActsExamples::EventAction::instance() {
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::EventAction::EventAction(std::reference_wrapper<std::vector<std::string>> processFilter) : G4UserEventAction(), m_processFilter(std::move(processFilter)) {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::EventAction::~EventAction() {
  s_instance = nullptr;
}

void ActsExamples::EventAction::BeginOfEventAction(const G4Event*) {
  SteppingAction::instance()->clear();
  m_event = HepMC3::GenEvent(HepMC3::Units::GEV, HepMC3::Units::MM);
}

void ActsExamples::EventAction::EndOfEventAction(const G4Event*) {
	// Walk over all vertices
	for(const auto& vertex : m_event->vertices())
	{
		// Consider only 1->1 vertices to maintain a correct history
		if(vertex->particles_in() == 1 && vertex->particles_out() == 1)
		{
			// Test for all attributes if one matches the filter pattern
			const std::vector<std::string> vertexAttributes = vertex->attribute_names();
			for(const auto& att : vertexAttributes)
			{
				const std::string process = vertex->attribute_as_string(att);
				if(std::find(m_processFilter.begin(), m_processFilter.end(), process) != m_processFilter.end())
				{
					auto& particleIn = vertex->particles_in()[0];
					auto prodVertex = particleIn->production_vertex();
					const int trackId = particleIn->attribute<HepMC3::IntAttribute>("TrackID");
					
					
					auto& particleOut = vertex->particles_out()[0];
					
					break;
				}
			}
		}
	}
}

void ActsExamples::EventAction::clear() {
  SteppingAction::instance()->clear();
}

HepMC3::GenEvent& ActsExamples::EventAction::event() {
  return m_event;
}