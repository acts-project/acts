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
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include "edm4hep/SimTrackerHitCollection.h"

namespace ActsExamples {

EDM4hepSimHitReader::EDM4hepSimHitReader(
    const EDM4hepSimHitReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepSimHitReader", level)) {
  m_reader.openFile(m_cfg.inputPath);
  m_store.setReader(&m_reader);

  m_eventsRange = std::make_pair(0, m_reader.getEntries());

  m_collections = m_reader.getCollectionIDTable()->names();
}

std::string EDM4hepSimHitReader::EDM4hepSimHitReader::name() const {
  return "EDM4hepSimHitReader";
}

std::pair<size_t, size_t> EDM4hepSimHitReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode EDM4hepSimHitReader::read(const AlgorithmContext& ctx) {
  SimHitContainer::sequence_type unordered;

  m_store.clear();
  m_reader.goToEvent(ctx.eventNumber);

  // TODO does it make sense to query all of them?
  for (const auto& name : m_collections) {
    auto& collection = m_store.get<podio::CollectionBase>(name);

    if (!collection.isValid()) {
      continue;
    }

    if (collection.getValueTypeName() == "edm4hep::SimTrackerHit") {
      for (const auto& simTrackerHit :
           (const edm4hep::SimTrackerHitCollection&)collection) {
        try {
          auto hit = EDM4hepUtil::fromSimHit(
              simTrackerHit,
              [](edm4hep::MCParticle particle) {
                ActsFatras::Barcode result;
                // TODO dont use podio internal id
                result.setParticle(particle.id());
                return result;
              },
              [&](std::uint64_t cellId) {
                auto detElement = m_cfg.dd4hepGeometryService->lcdd()
                                      ->volumeManager()
                                      .lookupDetElement(cellId);
                Acts::GeometryIdentifier result = detElement.volumeID();
                return result;
              });
          unordered.push_back(std::move(hit));
        } catch (...) {
          m_logger->log(Acts::Logging::Level::ERROR,
                        "EDM4hepSimHitReader: failed to convert SimTrackerHit");
          continue;
        }
      }
    }
  }

  m_reader.endOfEvent();

  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
