// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/ParametricProcessGenerator.hpp"

#include <random>

#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

FW::ParametricProcessGenerator::ParametricProcessGenerator(
    const FW::ParametricProcessGenerator::Config& cfg)
    : m_cfg(cfg),
      m_charge(ActsFatras::findCharge(m_cfg.pdg)),
      m_mass(ActsFatras::findMass(m_cfg.pdg)) {}

std::vector<FW::SimVertex> FW::ParametricProcessGenerator::operator()(
    FW::RandomEngine& rng) const {
  using UniformReal = std::uniform_real_distribution<double>;
  using UniformIndex = std::uniform_int_distribution<size_t>;

  UniformReal d0Dist(m_cfg.d0Range[0], m_cfg.d0Range[1]);
  UniformReal z0Dist(m_cfg.z0Range[0], m_cfg.z0Range[1]);
  UniformReal t0Dist(m_cfg.t0Range[0], m_cfg.t0Range[1]);
  UniformReal phiDist(m_cfg.phiRange[0], m_cfg.phiRange[1]);
  UniformReal etaDist(m_cfg.etaRange[0], m_cfg.etaRange[1]);
  UniformReal ptDist(m_cfg.ptRange[0], m_cfg.ptRange[1]);
  // choose between particle/anti-particle if requested
  UniformIndex particleTypeChoice(0u, m_cfg.randomizeCharge ? 1u : -0u);
  // (anti-)particle choice is one random draw but defines two properties
  const Acts::PdgParticle pdgChoices[] = {
      m_cfg.pdg,
      static_cast<Acts::PdgParticle>(-m_cfg.pdg),
  };
  const double qChoices[] = {m_charge, -m_charge};

  // create empty process vertex
  SimVertex process(SimVertex::Vector4::Zero());

  // counter will be reused as barcode particle number which must be non-zero.
  for (size_t ip = 1; ip <= m_cfg.numParticles; ++ip) {
    const auto d0 = d0Dist(rng);
    const auto z0 = z0Dist(rng);
    const auto t0 = t0Dist(rng);
    const auto phi = phiDist(rng);
    const auto eta = etaDist(rng);
    const auto pt = ptDist(rng);
    const auto type = particleTypeChoice(rng);
    const auto pdg = pdgChoices[type];
    const auto q = qChoices[type];

    // all particles are treated as originating from the same primary vertex
    const auto pid = ActsFatras::Barcode(0u).setParticle(ip);
    // construct the particle;
    ActsFatras::Particle particle(pid, pdg, q, m_mass);
    particle.setPosition4(d0 * std::sin(phi), d0 * -std::cos(phi), z0, t0);
    particle.setDirection(Acts::makeDirectionUnitFromPhiEta(phi, eta));
    // TODO check abs p value
    particle.setAbsMomentum(pt / std::cosh(eta));

    process.outgoing.push_back(std::move(particle));
  }

  return {process};
}
