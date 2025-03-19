// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitOutputConverter.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepSimHitOutputConverter::EDM4hepSimHitOutputConverter(
    const EDM4hepSimHitOutputConverter::Config& config,
    Acts::Logging::Level level)
    : EDM4hepOutputConverter("EDM4hepSimHitOutputConverter", level),
      m_cfg(config) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }

  if (m_cfg.outputSimTrackerHits.empty()) {
    throw std::invalid_argument("Missing output sim tracker hit name");
  }

  if (m_cfg.outputParticles.empty() != m_cfg.inputParticles.empty()) {
    throw std::invalid_argument(
        "Output particles and input particles must both be set or not set");
  }

  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_outputParticles.maybeInitialize(m_cfg.outputParticles);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_outputSimTrackerHits.initialize(m_cfg.outputSimTrackerHits);
}

ProcessCode EDM4hepSimHitOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  EDM4hepUtil::MapParticleIdTo particleMapper;
  std::unordered_map<ActsFatras::Barcode, edm4hep::MutableMCParticle>
      particleMap;

  edm4hep::SimTrackerHitCollection simTrackerHitCollection;

  if (!m_cfg.inputParticles.empty()) {
    edm4hep::MCParticleCollection mcParticles;
    auto particles = m_inputParticles(ctx);

    for (const auto& particle : particles) {
      auto p = mcParticles->create();
      particleMap[particle.particleId()] = p;
      EDM4hepUtil::writeParticle(particle, p);
    }

    particleMapper = [&](ActsFatras::Barcode particleId) {
      auto it = particleMap.find(particleId);
      if (it == particleMap.end()) {
        throw std::runtime_error("Particle not found in map");
      }
      return it->second;
    };
    m_outputParticles(ctx, std::move(mcParticles));
  }

  const auto& simHits = m_inputSimHits(ctx);

  for (const auto& simHit : simHits) {
    auto simTrackerHit = simTrackerHitCollection->create();
    EDM4hepUtil::writeSimHit(
        simHit, simTrackerHit, particleMapper,
        [](Acts::GeometryIdentifier id) { return id.value(); });
  }

  m_outputSimTrackerHits(ctx, std::move(simTrackerHitCollection));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepSimHitOutputConverter::collections() const {
  std::vector<std::string> result{m_cfg.outputSimTrackerHits};
  if (!m_cfg.outputParticles.empty()) {
    result.push_back(m_cfg.outputParticles);
  }
  return result;
}

}  // namespace ActsExamples
