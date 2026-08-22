// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/PrimaryGenerator.hpp"

#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include <cmath>
#include <numbers>

namespace ActsFatras::Synthetic::detail {

PrimaryGenerator::PrimaryGenerator(const EventConfig& config)
    : m_cfg(config), m_mass(config.particleMass()) {}

void PrimaryGenerator::generate(
    std::mt19937& rng, std::vector<Helix>& tracks,
    std::vector<GeneratedParticle>& particles) const {
  constexpr float pi = std::numbers::pi_v<float>;

  const GenerationConfig& cfg = m_cfg.generation;
  const std::size_t numPrimaries = cfg.numPrimaries();
  tracks.reserve(tracks.size() + numPrimaries);
  particles.reserve(particles.size() + numPrimaries);

  for (std::size_t t = 0; t < numPrimaries; ++t) {
    const float pt =
        samplePt(sampleUniform01(rng), cfg.minPt, cfg.ptScale, cfg.ptExponent);
    const float rapidity =
        sampleRapidity(rng, cfg.minRapidity, cfg.maxRapidity, cfg.rapidityEdge,
                       cfg.rapidityEdgeWidth);
    // cot(theta) = pz / pT = (mT / pT) sinh(y); the mass is the whole of the
    // difference between a rapidity and a pseudorapidity
    const float cotTheta = std::hypot(1.f, m_mass / pt) * std::sinh(rapidity);
    const float d0 = cfg.d0Sigma * sampleNormal(rng);
    const float phi = sampleUniform(rng, -pi, pi);
    const float z0 = cfg.beamspotSigmaZ * sampleNormal(rng);
    // as many of one sign as of the other, there being no charge asymmetry
    // worth describing in a minimum-bias event
    const float charge = sampleCharge(rng);
    const Helix helix =
        makeHelix(d0, phi, z0, cotTheta, pt, charge, m_cfg.bFieldZ);
    tracks.push_back(helix);

    // out of the luminous region, so generation zero and no production vertex
    // to speak of
    particles.push_back(makeParticle(helix, pt));
  }
}

}  // namespace ActsFatras::Synthetic::detail
