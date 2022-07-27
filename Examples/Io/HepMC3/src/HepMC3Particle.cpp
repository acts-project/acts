// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Particle.hpp"

#include "ActsExamples/Io/HepMC3/HepMC3Vertex.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

ActsExamples::SimParticle ActsExamples::HepMC3Particle::particle(
    HepMC3::ConstGenParticlePtr particle) {
  // TODO this is probably not quite right
  ActsFatras::Barcode particleId;
  particleId.setParticle(particle->id());
  Acts::PdgParticle pdg = static_cast<Acts::PdgParticle>(particle->pid());
  SimParticle fw(particleId, static_cast<Acts::PdgParticle>(particle->pid()),
                 ActsFatras::findCharge(pdg), particle->generated_mass());
  fw.setDirection(particle->momentum().x(), particle->momentum().y(),
                  particle->momentum().z());
  fw.setAbsoluteMomentum(particle->momentum().p3mod());
  return fw;
}

int ActsExamples::HepMC3Particle::id(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  return particle->id();
}

std::unique_ptr<ActsExamples::SimVertex>
ActsExamples::HepMC3Particle::productionVertex(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  // Return the vertex if it exists
  if (particle->production_vertex()) {
    return HepMC3Vertex::processVertex(
        std::make_shared<HepMC3::GenVertex>(*particle->production_vertex()));
  } else {
    return nullptr;
  }
}

std::unique_ptr<ActsExamples::SimVertex>
ActsExamples::HepMC3Particle::endVertex(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  // Return the vertex if it exists
  if (particle->end_vertex()) {
    return HepMC3Vertex::processVertex(
        std::make_shared<HepMC3::GenVertex>(*(particle->end_vertex())));
  } else {
    return nullptr;
  }
}

int ActsExamples::HepMC3Particle::pdgID(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  return particle->pid();
}

Acts::Vector3 ActsExamples::HepMC3Particle::momentum(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  Acts::Vector3 mom;
  mom(0) = particle->momentum().x();
  mom(1) = particle->momentum().y();
  mom(2) = particle->momentum().z();
  return mom;
}

double ActsExamples::HepMC3Particle::energy(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  return particle->momentum().e();
}

double ActsExamples::HepMC3Particle::mass(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  return particle->generated_mass();
}

double ActsExamples::HepMC3Particle::charge(
    const std::shared_ptr<HepMC3::GenParticle> particle) {
  return ActsFatras::findCharge(
      static_cast<Acts::PdgParticle>(particle->pid()));
}

void ActsExamples::HepMC3Particle::pdgID(
    std::shared_ptr<HepMC3::GenParticle> particle, const int pid) {
  particle->set_pid(pid);
}

void ActsExamples::HepMC3Particle::momentum(
    std::shared_ptr<HepMC3::GenParticle> particle, const Acts::Vector3& mom) {
  HepMC3::FourVector fVec(mom(0), mom(1), mom(2), particle->momentum().e());
  particle->set_momentum(fVec);
}

void ActsExamples::HepMC3Particle::energy(
    std::shared_ptr<HepMC3::GenParticle> particle, const double energy) {
  HepMC3::FourVector fVec(particle->momentum().x(), particle->momentum().y(),
                          particle->momentum().z(), energy);
  particle->set_momentum(fVec);
}

void ActsExamples::HepMC3Particle::mass(
    std::shared_ptr<HepMC3::GenParticle> particle, const double mass) {
  particle->set_generated_mass(mass);
}
