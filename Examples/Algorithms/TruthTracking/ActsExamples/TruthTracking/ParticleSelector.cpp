// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSelector.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
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
  if (!m_cfg.outputParticlesFinal.empty() &&
      m_cfg.inputParticlesFinal.empty()) {
    throw std::invalid_argument(
        "Output final particles collection requires input final particles");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticlesFinal.maybeInitialize(m_cfg.inputParticlesFinal);
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputParticlesFinal.maybeInitialize(m_cfg.outputParticlesFinal);

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
  // prepare input/ output types
  const SimParticleContainer& inputParticles = m_inputParticles(ctx);
  const SimParticleContainer& inputParticlesFinal =
      (m_inputParticlesFinal.isInitialized()) ? m_inputParticlesFinal(ctx)
                                              : inputParticles;

  std::size_t nInvalidCharge = 0;
  std::size_t nInvalidMeasurementCount = 0;

  // helper functions to select tracks
  auto within = [](auto x, auto min, auto max) {
    return (min <= x) && (x < max);
  };

  auto isValidParticle = [&](const ActsFatras::Particle& p) {
    const auto eta = Acts::VectorHelpers::eta(p.direction());
    const auto phi = Acts::VectorHelpers::phi(p.direction());
    const auto rho = Acts::VectorHelpers::perp(p.position());
    // define charge selection
    const bool validNeutral = (p.charge() == 0) && !m_cfg.removeNeutral;
    const bool validCharged = (p.charge() != 0) && !m_cfg.removeCharged;
    const bool validCharge = validNeutral || validCharged;
    const bool validSecondary = !m_cfg.removeSecondaries || !p.isSecondary();

    nInvalidCharge += static_cast<std::size_t>(!validCharge);

    bool validMeasurementCount = true;
    if (auto finalParticleIt = inputParticlesFinal.find(p.particleId());
        finalParticleIt != inputParticlesFinal.end()) {
      validMeasurementCount =
          within(finalParticleIt->numberOfHits(), m_cfg.measurementsMin,
                 m_cfg.measurementsMax);
    } else {
      ACTS_WARNING("No final particle found for " << p.particleId());
    }

    nInvalidMeasurementCount +=
        static_cast<std::size_t>(!validMeasurementCount);

    // Pdg selection
    bool validPdg = true;
    for (auto pdg : m_cfg.excludeAbsPdgs) {
      if (p.absolutePdg() == std::abs(pdg)) {
        validPdg = false;
        break;
      }
    }

    return validPdg && validCharge && validSecondary && validMeasurementCount &&
           within(p.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) &&
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) &&
           within(eta, m_cfg.etaMin, m_cfg.etaMax) &&
           within(phi, m_cfg.phiMin, m_cfg.phiMax) &&
           within(std::abs(p.position()[Acts::ePos2]), m_cfg.absZMin,
                  m_cfg.absZMax) &&
           within(rho, m_cfg.rhoMin, m_cfg.rhoMax) &&
           within(p.time(), m_cfg.timeMin, m_cfg.timeMax) &&
           within(p.mass(), m_cfg.mMin, m_cfg.mMax);
  };

  SimParticleContainer outputParticles;
  outputParticles.reserve(inputParticles.size());

  SimParticleContainer outputParticlesFinal;
  if (m_outputParticlesFinal.isInitialized()) {
    outputParticlesFinal.reserve(inputParticles.size());
  }

  // copy selected particles
  for (const auto& inputParticle : inputParticles) {
    if (!isValidParticle(inputParticle)) {
      continue;
    }

    outputParticles.insert(outputParticles.end(), inputParticle);

    if (m_outputParticlesFinal.isInitialized()) {
      if (auto particleFinalIt =
              inputParticlesFinal.find(inputParticle.particleId());
          particleFinalIt != inputParticlesFinal.end()) {
        outputParticlesFinal.insert(outputParticlesFinal.end(),
                                    *particleFinalIt);
      }
    }
  }
  outputParticles.shrink_to_fit();
  outputParticlesFinal.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputParticles.size() << " from "
                      << inputParticles.size() << " particles");
  ACTS_DEBUG("filtered out because of charge: " << nInvalidCharge);
  ACTS_DEBUG("filtered out because of measurement count: "
             << nInvalidMeasurementCount);

  m_outputParticles(ctx, std::move(outputParticles));
  if (m_outputParticlesFinal.isInitialized()) {
    m_outputParticlesFinal(ctx, std::move(outputParticlesFinal));
  }

  return ProcessCode::SUCCESS;
}
