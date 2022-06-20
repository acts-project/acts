// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHit.h"

namespace ActsExamples {

EDM4hepSimHitWriter::EDM4hepSimHitWriter(
    const EDM4hepSimHitWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputSimHits, "CsvSimHitWriter", level),
      m_cfg(config),
      m_writer(config.outputPath, &m_store) {
  ACTS_VERBOSE("Created output file " << config.outputPath);

  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }

  m_mcParticleCollection =
      &m_store.create<edm4hep::MCParticleCollection>(m_cfg.outputParticles);
  m_writer.registerForWrite(m_cfg.outputParticles);

  m_simTrackerHitCollection = &m_store.create<edm4hep::SimTrackerHitCollection>(
      m_cfg.outputSimTrackerHits);
  m_writer.registerForWrite(m_cfg.outputSimTrackerHits);
}

ActsExamples::ProcessCode EDM4hepSimHitWriter::endRun() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepSimHitWriter::writeT(const AlgorithmContext& ctx,
                                        const SimHitContainer& simHits) {
  EDM4hepUtil::MapParticleIdTo particleMapper;
  std::unordered_map<ActsFatras::Barcode, edm4hep::MutableMCParticle>
      particleMap;

  if (!m_cfg.inputParticles.empty()) {
    auto particles =
        ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

    for (const auto& particle : particles) {
      auto p = m_mcParticleCollection->create();
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
  }

  for (const auto& simHit : simHits) {
    auto simTrackerHit = m_simTrackerHitCollection->create();
    EDM4hepUtil::writeSimHit(
        simHit, simTrackerHit, particleMapper,
        [](Acts::GeometryIdentifier id) { return id.value(); });
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
