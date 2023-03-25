// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepParticleReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

namespace ActsExamples {

EDM4hepParticleReader::EDM4hepParticleReader(
    const EDM4hepParticleReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepParticleReader", level)) {
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_reader.openFile(m_cfg.inputPath);
  m_store.setReader(&m_reader);

  m_eventsRange = std::make_pair(0, m_reader.getEntries());

  m_mcParticleCollection =
      &m_store.get<edm4hep::MCParticleCollection>(m_cfg.inputParticles);

  m_outputParticles.initialize(m_cfg.outputParticles);
}

std::string EDM4hepParticleReader::name() const {
  return "EDM4hepParticleReader";
}

std::pair<size_t, size_t> EDM4hepParticleReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode EDM4hepParticleReader::read(const AlgorithmContext& ctx) {
  m_store.clear();
  m_reader.goToEvent(ctx.eventNumber);

  SimParticleContainer::sequence_type unordered;

  for (const auto& mcParticle : *m_mcParticleCollection) {
    auto particle =
        EDM4hepUtil::readParticle(mcParticle, [](const edm4hep::MCParticle& p) {
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

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
