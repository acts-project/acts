// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSelector.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>

ActsExamples::ParticleSelector::ParticleSelector(const Config& config,
                                                 Acts::Logging::Level level)
    : IAlgorithm("ParticleSelector", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputParticles.initialize(m_cfg.outputParticles);

  ACTS_DEBUG("selection particle rho [" << m_cfg.rhoMin << "," << m_cfg.rhoMax
                                        << ")");
  ACTS_DEBUG("selection particle |z| [" << m_cfg.absZMin << "," << m_cfg.absZMax
                                        << ")");
  ACTS_DEBUG("selection particle time [" << m_cfg.timeMin << ","
                                         << m_cfg.timeMax << ")");
  ACTS_DEBUG("selection particle phi [" << m_cfg.phiMin << "," << m_cfg.phiMax
                                        << ")");
  ACTS_DEBUG("selection particle eta [" << m_cfg.etaMin << "," << m_cfg.etaMax
                                        << ")");
  ACTS_DEBUG("selection particle |eta| [" << m_cfg.absEtaMin << ","
                                          << m_cfg.absEtaMax << ")");
  ACTS_DEBUG("selection particle pt [" << m_cfg.ptMin << "," << m_cfg.ptMax
                                       << ")");
  ACTS_DEBUG("selection particle m [" << m_cfg.mMin << "," << m_cfg.mMax
                                      << ")");
  ACTS_DEBUG("remove charged particles " << m_cfg.removeCharged);
  ACTS_DEBUG("remove neutral particles " << m_cfg.removeNeutral);
  ACTS_DEBUG("remove secondary particles " << m_cfg.removeSecondaries);
}

ActsExamples::ProcessCode ActsExamples::ParticleSelector::execute(
    const AlgorithmContext& ctx) const {
  // helper functions to select tracks
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidParticle = [&](const ActsFatras::Particle& p) {
    const auto eta = Acts::VectorHelpers::eta(p.direction());
    const auto phi = Acts::VectorHelpers::phi(p.direction());
    const auto rho = Acts::VectorHelpers::perp(p.position());
    // define charge selection
    const bool validNeutral = (p.charge() == 0) and not m_cfg.removeNeutral;
    const bool validCharged = (p.charge() != 0) and not m_cfg.removeCharged;
    const bool validCharge = validNeutral or validCharged;
    const bool validSecondary = not m_cfg.removeSecondaries or !p.isSecondary();
    return validCharge and validSecondary and
           within(p.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(phi, m_cfg.phiMin, m_cfg.phiMax) and
           within(std::abs(p.position()[Acts::ePos2]), m_cfg.absZMin,
                  m_cfg.absZMax) and
           within(rho, m_cfg.rhoMin, m_cfg.rhoMax) and
           within(p.time(), m_cfg.timeMin, m_cfg.timeMax) and
           within(p.mass(), m_cfg.mMin, m_cfg.mMax);
  };

  // prepare input/ output types
  const auto& inputParticles = m_inputParticles(ctx);

  SimParticleContainer outputParticles;
  outputParticles.reserve(inputParticles.size());

  // copy selected particles
  for (const auto& inputParticle : inputParticles) {
    if (isValidParticle(inputParticle)) {
      // the input parameters should already be
      outputParticles.insert(outputParticles.end(), inputParticle);
    }
  }
  outputParticles.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputParticles.size() << " from "
                      << inputParticles.size() << " particles");

  m_outputParticles(ctx, std::move(outputParticles));
  return ProcessCode::SUCCESS;
}
