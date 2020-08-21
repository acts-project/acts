// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Plugins/HepMC3/HepMC3Vertex.hpp"

#include "ActsExamples/Plugins/HepMC3/HepMC3Particle.hpp"

std::vector<ActsExamples::SimParticle>
ActsExamples::HepMC3Vertex::genParticlesToActs(
    const std::vector<HepMC3::GenParticlePtr>& genParticles) {
  HepMC3Particle simPart;

  std::vector<SimParticle> actsParticles;
  // Translate all particles
  for (auto& genParticle : genParticles)
    actsParticles.push_back(*(
        simPart.particle(std::make_shared<HepMC3::GenParticle>(*genParticle))));
  return actsParticles;
}

std::unique_ptr<ActsExamples::SimVertex>
ActsExamples::HepMC3Vertex::processVertex(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  SimVertex vtx({vertex->position().x(), vertex->position().y(),
                 vertex->position().z(), vertex->position().t()});
  vtx.incoming = genParticlesToActs(vertex->particles_in());
  vtx.outgoing = genParticlesToActs(vertex->particles_out());
  // Create Acts vertex
  return std::make_unique<SimVertex>(std::move(vtx));
}

bool ActsExamples::HepMC3Vertex::inEvent(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  return vertex->in_event();
}

int ActsExamples::HepMC3Vertex::id(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  return vertex->id();
}

std::vector<ActsExamples::SimParticle> ActsExamples::HepMC3Vertex::particlesIn(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  return genParticlesToActs(vertex->particles_in());
}

std::vector<ActsExamples::SimParticle> ActsExamples::HepMC3Vertex::particlesOut(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  return genParticlesToActs(vertex->particles_out());
}

Acts::Vector3D ActsExamples::HepMC3Vertex::position(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  Acts::Vector3D vec;
  vec(0) = vertex->position().x();
  vec(1) = vertex->position().y();
  vec(2) = vertex->position().z();
  return vec;
}

double ActsExamples::HepMC3Vertex::time(
    const std::shared_ptr<HepMC3::GenVertex> vertex) {
  return vertex->position().t();
}

HepMC3::GenParticlePtr ActsExamples::HepMC3Vertex::actsParticleToGen(
    std::shared_ptr<SimParticle> actsParticle) {
  // Extract momentum and energy from Acts particle for HepMC3::FourVector
  const auto mom = actsParticle->momentum4();
  const HepMC3::FourVector vec(mom[0], mom[1], mom[2], mom[3]);
  // Create HepMC3::GenParticle
  HepMC3::GenParticle genParticle(vec, actsParticle->pdg());
  genParticle.set_generated_mass(actsParticle->mass());

  return std::shared_ptr<HepMC3::GenParticle>(&genParticle);
}

void ActsExamples::HepMC3Vertex::addParticleIn(
    std::shared_ptr<HepMC3::GenVertex> vertex,
    std::shared_ptr<SimParticle> particle) {
  vertex->add_particle_in(actsParticleToGen(particle));
}

void ActsExamples::HepMC3Vertex::addParticleOut(
    std::shared_ptr<HepMC3::GenVertex> vertex,
    std::shared_ptr<SimParticle> particle) {
  vertex->add_particle_out(actsParticleToGen(particle));
}

HepMC3::GenParticlePtr ActsExamples::HepMC3Vertex::matchParticles(
    const std::vector<HepMC3::GenParticlePtr>& genParticles,
    std::shared_ptr<SimParticle> actsParticle) {
  const auto id = actsParticle->particleId();
  // Search HepMC3::GenParticle with the same id as the Acts particle
  for (auto& genParticle : genParticles) {
    if (genParticle->id() == id) {
      // Return particle if found
      return genParticle;
    }
  }
  return nullptr;
}

void ActsExamples::HepMC3Vertex::removeParticleIn(
    std::shared_ptr<HepMC3::GenVertex> vertex,
    std::shared_ptr<SimParticle> particle) {
  // Remove particle if it exists
  if (HepMC3::GenParticlePtr genParticle =
          matchParticles(vertex->particles_in(), particle))
    vertex->remove_particle_in(genParticle);
}

void ActsExamples::HepMC3Vertex::removeParticleOut(
    std::shared_ptr<HepMC3::GenVertex> vertex,
    std::shared_ptr<SimParticle> particle) {
  // Remove particle if it exists
  if (HepMC3::GenParticlePtr genParticle =
          matchParticles(vertex->particles_out(), particle))
    vertex->remove_particle_out(genParticle);
}

void ActsExamples::HepMC3Vertex::position(
    const std::shared_ptr<HepMC3::GenVertex> vertex, Acts::Vector3D pos) {
  HepMC3::FourVector fVec(pos(0), pos(1), pos(2), vertex->position().t());
  vertex->set_position(fVec);
}

void ActsExamples::HepMC3Vertex::time(
    const std::shared_ptr<HepMC3::GenVertex> vertex, double time) {
  HepMC3::FourVector fVec(vertex->position().x(), vertex->position().y(),
                          vertex->position().z(), time);
  vertex->set_position(fVec);
}
