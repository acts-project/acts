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

namespace {

bool findAttribute(HepMC3::ConstGenVertexPtr vertex, const std::vector<std::string>& processFilter) {
	
	// Consider only 1->1 vertices to maintain a correct history
	if(vertex->particles_in().size() == 1 && vertex->particles_out().size() == 1)
	{
for(const auto& s : vertex->attribute_names())
	std::cout << s << " ";
std::cout << std::endl;	
		// Test for all attributes if one matches the filter pattern
		const std::vector<std::string> vertexAttributes = vertex->attribute_names();
		for(const auto& att : vertexAttributes)
		{
			const std::string process = vertex->attribute_as_string(att);
			if(std::find(processFilter.begin(), processFilter.end(), process) != processFilter.end())
			{
				return true;
			}
		}
	}
	return false;
}

void reduceVertex(HepMC3::GenEvent& event, HepMC3::GenVertexPtr vertex, const std::vector<std::string>& processFilter) {
	auto particleIn = vertex->particles_in()[0];
	auto particleOut = vertex->particles_out()[0];
	event.remove_vertex(vertex);

	auto reducedVertex = std::make_shared<HepMC3::GenVertex>();
	event.add_vertex(reducedVertex);
	
	while(findAttribute(particleIn->production_vertex(), processFilter))
	{
		auto nextParticle = particleIn->production_vertex()->particles_in()[0];
		event.remove_particle(particleIn);
		particleIn = nextParticle;
		event.remove_vertex(particleIn->end_vertex());
	}
	while(findAttribute(particleOut->end_vertex(), processFilter))
	{	
		auto nextParticle = particleOut->end_vertex()->particles_out()[0];
		event.remove_particle(particleOut);
		particleOut = nextParticle;
		event.remove_vertex(particleOut->end_vertex());		
	}		
	reducedVertex->add_particle_in(particleIn);
	reducedVertex->add_particle_out(particleOut);
	vertex = reducedVertex;
}

void followOutgoingParticles(HepMC3::GenEvent& event, HepMC3::GenVertexPtr vertex, const std::vector<std::string>& processFilter)
{
	if(findAttribute(vertex, processFilter))
	{
		reduceVertex(event, vertex, processFilter);
	}
	for(const auto& particle : vertex->particles_out())
	{
		followOutgoingParticles(event, particle->end_vertex(), processFilter);
	}
}
}

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
	
	if(m_event.vertices().empty())
		return;
		
	auto currentVertex = m_event.vertices()[0];
	followOutgoingParticles(m_event, currentVertex, m_processFilter);
}

void ActsExamples::EventAction::clear() {
  SteppingAction::instance()->clear();
}

HepMC3::GenEvent& ActsExamples::EventAction::event() {
  return m_event;
}