// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/HepMC3/HepMC3Event.hpp"

#include "ACTFW/Plugins/HepMC3/HepMC3Particle.hpp"
#include "ACTFW/Plugins/HepMC3/HepMC3Vertex.hpp"

///
/// Setter
///

void FW::HepMC3Event::momentumUnit(std::shared_ptr<HepMC3::GenEvent> event,
                                   const double momentumUnit) {
  // Check, if the momentum unit fits Acts::units::_MeV or _GeV
  HepMC3::Units::MomentumUnit mom;
  if (momentumUnit == Acts::units::_MeV)
    mom = HepMC3::Units::MomentumUnit::MEV;
  else if (momentumUnit == Acts::units::_GeV)
    mom = HepMC3::Units::MomentumUnit::GEV;
  else {
    // Report invalid momentum unit and set GeV
    std::cout << "Invalid unit of momentum: " << momentumUnit << std::endl;
    std::cout << "Momentum unit [GeV] will be used instead" << std::endl;
    mom = HepMC3::Units::MomentumUnit::GEV;
  }
  // Set units
  event->set_units(mom, event->length_unit());
}

void FW::HepMC3Event::lengthUnit(std::shared_ptr<HepMC3::GenEvent> event,
                                 const double lengthUnit) {
  // Check, if the length unit fits Acts::units::_mm or _cm
  HepMC3::Units::LengthUnit len;
  if (lengthUnit == Acts::units::_mm)
    len = HepMC3::Units::LengthUnit::MM;
  else if (lengthUnit == Acts::units::_cm)
    len = HepMC3::Units::LengthUnit::CM;
  else {
    // Report invalid length unit and set mm
    std::cout << "Invalid unit of length: " << lengthUnit << std::endl;
    std::cout << "Length unit [mm] will be used instead" << std::endl;
    len = HepMC3::Units::LengthUnit::MM;
  }

  // Set units
  event->set_units(event->momentum_unit(), len);
}

void FW::HepMC3Event::shiftPositionBy(std::shared_ptr<HepMC3::GenEvent> event,
                                      const Acts::Vector3D& deltaPos,
                                      const double deltaTime) {
  // Create HepMC3::FourVector from position and time for shift
  const HepMC3::FourVector vec(deltaPos(0), deltaPos(1), deltaPos(2),
                               deltaTime);
  event->shift_position_by(vec);
}

void FW::HepMC3Event::shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                                      const Acts::Vector3D& pos,
                                      const double time) {
  // Create HepMC3::FourVector from position and time for the new position
  const HepMC3::FourVector vec(pos(0), pos(1), pos(2), time);
  event->shift_position_to(vec);
}

void FW::HepMC3Event::shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                                      const Acts::Vector3D& pos) {
  // Create HepMC3::FourVector from position and time for the new position
  const HepMC3::FourVector vec(pos(0), pos(1), pos(2), event->event_pos().t());
  event->shift_position_to(vec);
}

void FW::HepMC3Event::shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                                      const double time) {
  // Create HepMC3::FourVector from position and time for the new position
  const HepMC3::FourVector vec(event->event_pos().x(), event->event_pos().y(),
                               event->event_pos().z(), time);
  event->shift_position_to(vec);
}

///
/// Adder
///

HepMC3::GenParticlePtr FW::HepMC3Event::actsParticleToGen(
    std::shared_ptr<SimParticle> actsParticle) {
  // Extract momentum and energy from Acts particle for HepMC3::FourVector
  const auto mom4 = actsParticle->momentum4();
  const HepMC3::FourVector vec(mom4[0], mom4[1], mom4[2], mom4[3]);
  // Create HepMC3::GenParticle
  HepMC3::GenParticle genParticle(vec, actsParticle->pdg());
  genParticle.set_generated_mass(actsParticle->mass());

  return std::shared_ptr<HepMC3::GenParticle>(&genParticle);
}

void FW::HepMC3Event::addParticle(std::shared_ptr<HepMC3::GenEvent> event,
                                  std::shared_ptr<SimParticle> particle) {
  // Add new particle
  event->add_particle(actsParticleToGen(particle));
}

HepMC3::GenVertexPtr FW::HepMC3Event::createGenVertex(
    const std::shared_ptr<SimVertex>& actsVertex) {
  const HepMC3::FourVector vec(
      actsVertex->position4[0], actsVertex->position4[1],
      actsVertex->position4[2], actsVertex->position4[3]);

  // Create vertex
  HepMC3::GenVertex genVertex(vec);

  // Store incoming particles
  for (auto& particle : actsVertex->incoming) {
    HepMC3::GenParticlePtr genParticle =
        actsParticleToGen(std::make_shared<SimParticle>(particle));
    genVertex.add_particle_in(genParticle);
  }
  // Store outgoing particles
  for (auto& particle : actsVertex->outgoing) {
    HepMC3::GenParticlePtr genParticle =
        actsParticleToGen(std::make_shared<SimParticle>(particle));
    genVertex.add_particle_out(genParticle);
  }
  return std::shared_ptr<HepMC3::GenVertex>(&genVertex);
}

void FW::HepMC3Event::addVertex(std::shared_ptr<HepMC3::GenEvent> event,
                                const std::shared_ptr<SimVertex> vertex) {
  // Add new vertex
  event->add_vertex(createGenVertex(vertex));
}

///
/// Remover
///

void FW::HepMC3Event::removeParticle(
    std::shared_ptr<HepMC3::GenEvent> event,
    const std::shared_ptr<SimParticle>& particle) {
  const std::vector<HepMC3::GenParticlePtr> genParticles = event->particles();
  const auto id = particle->particleId();
  // Search HepMC3::GenParticle with the same id as the Acts particle
  for (auto& genParticle : genParticles) {
    if (genParticle->id() == id) {
      // Remove particle if found
      event->remove_particle(genParticle);
      break;
    }
  }
}

bool FW::HepMC3Event::compareVertices(
    const std::shared_ptr<SimVertex>& actsVertex,
    const HepMC3::GenVertexPtr& genVertex) {
  // Compare position, time, number of incoming and outgoing particles between
  // both vertices. Return false if one criterium does not match, else true.
  HepMC3::FourVector genVec = genVertex->position();
  if (actsVertex->position4[0] != genVec.x())
    return false;
  if (actsVertex->position4[1] != genVec.y())
    return false;
  if (actsVertex->position4[2] != genVec.z())
    return false;
  if (actsVertex->position4[3] != genVec.t())
    return false;
  if (actsVertex->incoming.size() != genVertex->particles_in().size())
    return false;
  if (actsVertex->outgoing.size() != genVertex->particles_out().size())
    return false;
  return true;
}

void FW::HepMC3Event::removeVertex(std::shared_ptr<HepMC3::GenEvent> event,
                                   const std::shared_ptr<SimVertex>& vertex) {
  const std::vector<HepMC3::GenVertexPtr> genVertices = event->vertices();
  // Walk over every recorded vertex
  for (auto& genVertex : genVertices)
    if (compareVertices(vertex, genVertex)) {
      // Remove vertex if it matches actsVertex
      event->remove_vertex(genVertex);
      break;
    }
}

///
/// Getter
///

double FW::HepMC3Event::momentumUnit(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  // HepMC allows only MEV and GEV. This allows an easy identification.
  return (event->momentum_unit() == HepMC3::Units::MomentumUnit::MEV
              ? Acts::units::_MeV
              : Acts::units::_GeV);
}

double FW::HepMC3Event::lengthUnit(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  // HepMC allows only MM and CM. This allows an easy identification.
  return (event->length_unit() == HepMC3::Units::LengthUnit::MM
              ? Acts::units::_mm
              : Acts::units::_cm);
}

Acts::Vector3D FW::HepMC3Event::eventPos(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  // Extract the position from HepMC3::FourVector
  Acts::Vector3D vec;
  vec(0) = event->event_pos().x();
  vec(1) = event->event_pos().y();
  vec(2) = event->event_pos().z();
  return vec;
}

double FW::HepMC3Event::eventTime(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  // Extract the time from HepMC3::FourVector
  return event->event_pos().t();
}

std::vector<std::unique_ptr<FW::SimParticle>> FW::HepMC3Event::particles(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  std::vector<std::unique_ptr<SimParticle>> actsParticles;
  const std::vector<HepMC3::GenParticlePtr> genParticles = event->particles();

  HepMC3Particle simPart;

  // Translate all particles
  for (auto& genParticle : genParticles)
    actsParticles.push_back(
        simPart.particle(std::make_shared<HepMC3::GenParticle>(*genParticle)));

  return actsParticles;
}

std::vector<std::unique_ptr<FW::SimVertex>> FW::HepMC3Event::vertices(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  std::vector<std::unique_ptr<SimVertex>> actsVertices;
  const std::vector<HepMC3::GenVertexPtr> genVertices = event->vertices();

  HepMC3Vertex simVert;

  // Translate all vertices
  for (auto& genVertex : genVertices) {
    actsVertices.push_back(
        simVert.processVertex(std::make_shared<HepMC3::GenVertex>(*genVertex)));
  }
  return actsVertices;
}

std::vector<std::unique_ptr<FW::SimParticle>> FW::HepMC3Event::beams(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  std::vector<std::unique_ptr<SimParticle>> actsBeams;
  const std::vector<HepMC3::GenParticlePtr> genBeams = event->beams();

  HepMC3Particle simPart;

  // Translate beam particles and store the result
  for (auto& genBeam : genBeams)
    actsBeams.push_back(
        simPart.particle(std::make_shared<HepMC3::GenParticle>(*genBeam)));
  return actsBeams;
}

std::vector<std::unique_ptr<FW::SimParticle>> FW::HepMC3Event::finalState(
    const std::shared_ptr<HepMC3::GenEvent> event) {
  std::vector<HepMC3::GenParticlePtr> particles = event->particles();
  std::vector<std::unique_ptr<SimParticle>> fState;

  HepMC3Particle simPart;

  // Walk over every vertex
  for (auto& particle : particles) {
    // Collect particles without end vertex
    if (!particle->end_vertex())
      fState.push_back(
          simPart.particle(std::make_shared<HepMC3::GenParticle>(*particle)));
  }
  return fState;
}
