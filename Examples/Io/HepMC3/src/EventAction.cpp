// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"
#include <stdexcept>
#include "SteppingAction.hpp"
#include <G4Event.hh>
#include <G4RunManager.hh>

ActsExamples::EventAction* ActsExamples::EventAction::s_instance = nullptr;

ActsExamples::EventAction*
ActsExamples::EventAction::instance()
{
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::EventAction::EventAction() : G4UserEventAction()
{
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::EventAction::~EventAction()
{
  s_instance = nullptr;
}

void
ActsExamples::EventAction::BeginOfEventAction(const G4Event*)
{
  SteppingAction::instance()->clear();
  m_event = std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
}

void
ActsExamples::EventAction::EndOfEventAction(const G4Event*)
{
	//~ std::cout << "Particles: " << m_event->particles().size() << " | " << "Vertices: " << m_event->vertices().size() << std::endl;
	//~ std::cout << "Number of steps: " << SteppingAction::instance()->counter() << std::endl;
	//~ SteppingAction::instance()->counter() = 0;
	
	//~ for(const auto& part : m_event->particles())
	//~ {	
		//~ if(!part->production_vertex())
			//~ std::cout << "Production vertex missing" << std::endl;
		//~ if(!part->end_vertex())
			//~ std::cout << "End vertex missing" << std::endl;
	//~ }
	//~ for(const auto& vert : m_event->vertices())
	//~ {
		//~ if(vert->particles_out().size() > 1)
		//~ {
			//~ if(!vert->particles_in().empty())
			//~ {
				//~ auto partIn = vert->particles_in()[0];
				//~ auto prodVertex = partIn->production_vertex();
				//~ for(const auto& s : prodVertex->attribute_names())
					//~ std::cout << s << " ";
				//~ std::cout << "Wanted: " << partIn->id() << std::endl;
				//~ if(!prodVertex->attribute_names().empty())
					//~ std::cout << prodVertex->attribute<HepMC3::StringAttribute>("NextProcessOf-" + std::to_string(partIn->id()))->value() << " leads to: ";
			//~ }
			//~ for(const auto& s : vert->attribute_names())
				//~ std::cout << typeid(vert->attribute<HepMC3::StringAttribute>(s)).name() << std::endl;

			//~ std::cout << ": " << vert->particles_in().size() << "(" 
				//~ << (vert->particles_in().size() > 0 ? vert->particles_in()[0]->pid() : 0) << ") -> " << vert->particles_out().size() << " ";
			//~ for(const auto& part : vert->particles_out())
				//~ std::cout << "(" << part->pid() << ") ";
				
			//~ std::cout << "| (" 	<< vert->position().x() << ", " << vert->position().y() << ", " << vert->position().z() << ")" << std::endl;
		//~ }
	//~ }
	//~ std::cout << "Checks done" << std::endl;
}

void
ActsExamples::EventAction::clear()
{
	m_event = nullptr;
	SteppingAction::instance()->clear();
}

std::shared_ptr<HepMC3::GenEvent>
ActsExamples::EventAction::event() const
{
  return m_event;
}