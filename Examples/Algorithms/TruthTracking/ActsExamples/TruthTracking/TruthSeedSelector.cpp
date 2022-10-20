// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

TruthSeedSelector::TruthSeedSelector(const Config& config,
                                     Acts::Logging::Level level)
    : BareAlgorithm("TruthSeedSelector", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output truth particles collection");
  }
}

ProcessCode TruthSeedSelector::execute(const AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  // prepare input collections
  const auto& inputParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  // prepare output collection
  SimParticleContainer selectedParticles;
  selectedParticles.reserve(inputParticles.size());

  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidparticle = [&](const auto& p) {
    const auto eta = Acts::VectorHelpers::eta(p.unitDirection());
    const auto phi = Acts::VectorHelpers::phi(p.unitDirection());
    const auto rho = Acts::VectorHelpers::perp(p.position());
    // find the corresponding hits for this particle
    const auto& hits = makeRange(particleHitsMap.equal_range(p.particleId()));
    // number of recorded hits
    size_t nHits = hits.size();
    return within(rho, 0., m_cfg.rhoMax) and
           within(p.position().z(), m_cfg.zMin, m_cfg.zMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(phi, m_cfg.phiMin, m_cfg.phiMax) and
           within(p.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
           within(nHits, m_cfg.nHitsMin, m_cfg.nHitsMax) and
           (m_cfg.keepNeutral or (p.charge() != 0));
  };

  // create prototracks for all input particles
  for (const auto& particle : inputParticles) {
    if (isValidparticle(particle)) {
      selectedParticles.insert(particle);
    }
  }

  ctx.eventStore.add(m_cfg.outputParticles, std::move(selectedParticles));
  return ProcessCode::SUCCESS;
}
