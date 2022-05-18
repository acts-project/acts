// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitReader.hpp"

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "ActsFatras/EventData/Hit.hpp"

#include "Acts/Definitions/Units.hpp"

#include "edm4hep/SimTrackerHit.h"
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

  m_simHitCollections = m_reader.getCollectionIDTable()->names();
}

std::string ActsExamples::EDM4hepSimHitReader::EDM4hepSimHitReader::name()
    const {
  return "EDM4hepSimHitReader";
}

std::pair<size_t, size_t> ActsExamples::EDM4hepSimHitReader::availableEvents()
    const {
  return m_eventsRange;
}

ActsFatras::Hit convertEDM4hepSimHit(
    const edm4hep::SimTrackerHit& sth,
    ActsExamples::DD4hep::DD4hepGeometryService &geometryService) {
  std::ofstream debug("/home/andreas/debug.txt", std::ios::app);

  debug << "EDM4hep convert cell id " << sth.getCellID() << std::endl;

  try {
    auto detElement = geometryService.lcdd()->volumeManager().lookupDetElement(sth.getCellID());
    debug << "EDM4hep found detElement " << detElement.volumeID() << std::endl;
  } catch (...) {
    debug << "EDM4hep detElement not found" << std::endl;
  }

  const auto geometryId = sth.getCellID(); // TODO
  ActsFatras::Barcode particleId;
  particleId.setParticle(sth.getMCParticle().id()); // TODO

  const auto mass = sth.getMCParticle().getMass();
  const Acts::ActsVector<3> momentum{
      sth.getMomentum().x * Acts::UnitConstants::GeV,
      sth.getMomentum().y * Acts::UnitConstants::GeV,
      sth.getMomentum().z * Acts::UnitConstants::GeV,
  };
  const auto energy = std::sqrt(momentum.squaredNorm() + mass * mass);

  ActsFatras::Hit::Vector4 pos4{
      sth.getPosition().x * Acts::UnitConstants::mm,
      sth.getPosition().y * Acts::UnitConstants::mm,
      sth.getPosition().z * Acts::UnitConstants::mm,
      sth.getTime() * Acts::UnitConstants::ns,
  };
  ActsFatras::Hit::Vector4 mom4{
      momentum.x(),
      momentum.y(),
      momentum.z(),
      energy,
  };
  ActsFatras::Hit::Vector4 delta4{
      0 * Acts::UnitConstants::GeV, 0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,  // sth.getEDep()
  };
  int32_t index = -1;

  return ActsFatras::Hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                         index);
}

ActsExamples::ProcessCode ActsExamples::EDM4hepSimHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  SimHitContainer::sequence_type unordered;

  m_store.clear();
  m_reader.goToEvent(ctx.eventNumber);

  for (const auto& name : m_simHitCollections) {
    auto& collection = m_store.get<podio::CollectionBase>(name);

    if (!collection.isValid()) {
      continue;
    }

    if (collection.getTypeName() == "edm4hep::SimTrackerHitCollection") {
      for (const auto& sth :
           (const edm4hep::SimTrackerHitCollection&)collection) {
        auto hit = convertEDM4hepSimHit(sth, *m_cfg.dd4hepGeometryService);
        unordered.push_back(std::move(hit));
      }
    }
  }

  m_reader.endOfEvent();

  // write the ordered data to the EventStore (according to geometry_id).
  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
