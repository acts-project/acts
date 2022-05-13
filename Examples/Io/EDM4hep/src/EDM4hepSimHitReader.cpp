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
#include "ActsExamples/Utilities/Paths.hpp"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

ActsExamples::EDM4hepSimHitReader::EDM4hepSimHitReader(
    const ActsExamples::EDM4hepSimHitReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepSimHitReader", level)) {
  std::cerr << "EDM4hepSimHitReader: input path " << m_cfg.inputPath
            << std::endl;

  m_reader.openFile(m_cfg.inputPath);
  m_store.setReader(&m_reader);

  m_eventsRange = std::make_pair(0, m_reader.getEntries());
}

std::string ActsExamples::EDM4hepSimHitReader::EDM4hepSimHitReader::name()
    const {
  return "EDM4hepSimHitReader";
}

std::pair<size_t, size_t> ActsExamples::EDM4hepSimHitReader::availableEvents()
    const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::EDM4hepSimHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  SimHitContainer::sequence_type unordered;

  m_store.clear();
  m_reader.goToEvent(ctx.eventNumber);

  for (const auto& name : m_reader.getCollectionIDTable()->names()) {
    auto& sths = m_store.get<edm4hep::SimTrackerHitCollection>(name);

    if (!sths.isValid()) {
      continue;
    }

    for (const auto& sth : sths) {
      const auto geometryId = 0; // Acts::GeometryIdentifier(sth.getCellID());  // TODO
      const auto particleId = 0; // ActsFatras::Barcode(sth.getMCParticle().getPDG());  // TODO

      ActsFatras::Hit::Vector4 pos4{
          sth.getPosition().x * Acts::UnitConstants::mm,
          sth.getPosition().y * Acts::UnitConstants::mm,
          sth.getPosition().z * Acts::UnitConstants::mm,
          sth.getTime() * Acts::UnitConstants::ns,
      };
      ActsFatras::Hit::Vector4 mom4{
          sth.getMomentum().x * Acts::UnitConstants::GeV,
          sth.getMomentum().y * Acts::UnitConstants::GeV,
          sth.getMomentum().z * Acts::UnitConstants::GeV,
          100 * Acts::UnitConstants::GeV,  // TODO
      };
      ActsFatras::Hit::Vector4 delta4{
          0 * Acts::UnitConstants::GeV,  // TODO
          0 * Acts::UnitConstants::GeV,  // TODO
          0 * Acts::UnitConstants::GeV,  // TODO
          sth.getEDep() * Acts::UnitConstants::GeV,
      };
      int32_t index = -1;

      ActsFatras::Hit hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                          index);
      unordered.push_back(std::move(hit));
    }
  }

  m_reader.endOfEvent();

  // write the ordered data to the EventStore (according to geometry_id).
  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
