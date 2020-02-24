// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/Particle.hpp"

#include "ActsFatras/Utilities/ParticleData.hpp"

ActsFatras::Particle::Particle(Barcode particleId, Acts::PdgParticle pdg)
    : Particle(particleId, pdg, findCharge(pdg), findMass(pdg)) {}

std::ostream& ActsFatras::operator<<(std::ostream& os,
                                     const ActsFatras::Particle& particle) {
  os << "particle_id=" << particle.particleId();
  os << " process=" << particle.process();
  os << " pdg=" << particle.pdg();
  os << " q=" << particle.charge();
  os << " m=" << particle.mass();
  return os;
}
