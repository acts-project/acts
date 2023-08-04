// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepSimHitReader::EDM4hepSimHitReader(
    const EDM4hepSimHitReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepSimHitReader", level)) {
  m_reader.openFile(m_cfg.inputPath);

  m_outputParticles.maybeInitialize(m_cfg.outputParticles);
  m_outputSimHits.initialize(m_cfg.outputSimHits);

  m_eventsRange = std::make_pair(0, m_reader.getEntries("events"));
}

std::string EDM4hepSimHitReader::EDM4hepSimHitReader::name() const {
  return "EDM4hepSimHitReader";
}

std::pair<size_t, size_t> EDM4hepSimHitReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode EDM4hepSimHitReader::read(const AlgorithmContext& ctx) {
  podio::Frame frame = m_reader.readEntry("events", ctx.eventNumber);

  const auto& mcParticleCollection =
      frame.get<edm4hep::MCParticleCollection>(m_cfg.inputParticles);

  if (!m_cfg.outputParticles.empty()) {
    SimParticleContainer::sequence_type unordered;

    for (const auto& mcParticle : mcParticleCollection) {
      auto particle = EDM4hepUtil::readParticle(
          mcParticle, [](const edm4hep::MCParticle& p) {
            ActsFatras::Barcode result;
            // TODO dont use podio internal id
            result.setParticle(p.id());
            return result;
          });
      unordered.push_back(particle);
    }

    // Write ordered particles container to the EventStore
    SimParticleContainer particles;
    particles.insert(unordered.begin(), unordered.end());
    m_outputParticles(ctx, std::move(particles));
  }

  SimHitContainer::sequence_type unordered;

  const auto& simTrackerHitCollection =
      frame.get<edm4hep::SimTrackerHitCollection>(m_cfg.inputSimTrackerHits);

  for (const auto& simTrackerHit : simTrackerHitCollection) {
    try {
      auto hit = EDM4hepUtil::readSimHit(
          simTrackerHit,
          [](const edm4hep::MCParticle& particle) {
            ActsFatras::Barcode result;
            // TODO dont use podio internal id
            result.setParticle(particle.id());
            return result;
          },
          [&](std::uint64_t cellId) {
            auto detElement = m_cfg.dd4hepDetector->geometryService->detector()
                                  .volumeManager()
                                  .lookupDetElement(cellId);
            Acts::GeometryIdentifier result = detElement.volumeID();
            return result;
          });
      unordered.push_back(std::move(hit));
    } catch (...) {
      ACTS_ERROR("EDM4hepSimHitReader: failed to convert SimTrackerHit");
      continue;
    }
  }

  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());
  m_outputSimHits(ctx, std::move(simHits));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
