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
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

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

  m_simTrackerHitCollection =
      &m_store.create<edm4hep::SimTrackerHitCollection>("ActsSimTrackerHits");
  m_writer.registerForWrite("ActsSimTrackerHits");
}

EDM4hepSimHitWriter::~EDM4hepSimHitWriter() {
  m_writer.finish();
}

ProcessCode EDM4hepSimHitWriter::writeT(const AlgorithmContext& ctx,
                                        const SimHitContainer& simHits) {
  for (const auto& simHit : simHits) {
    auto simTrackerHit = m_simTrackerHitCollection->create();

    // local simhit information in global coord.
    const Acts::Vector4& globalPos4 = simHit.fourPosition();
    const Acts::Vector4& momentum4Before = simHit.momentum4Before();
    const auto delta4 = simHit.momentum4After() - momentum4Before;

    // TODO set particle

    simTrackerHit.setCellID(simHit.geometryId().value());

    simTrackerHit.setTime(globalPos4[Acts::eTime] / Acts::UnitConstants::ns);
    simTrackerHit.setPosition({
        globalPos4[Acts::ePos0] / Acts::UnitConstants::mm,
        globalPos4[Acts::ePos1] / Acts::UnitConstants::mm,
        globalPos4[Acts::ePos2] / Acts::UnitConstants::mm,
    });

    simTrackerHit.setMomentum({
        (float)(momentum4Before[Acts::eMom0] / Acts::UnitConstants::GeV),
        (float)(momentum4Before[Acts::eMom1] / Acts::UnitConstants::GeV),
        (float)(momentum4Before[Acts::eMom2] / Acts::UnitConstants::GeV),
    });

    simTrackerHit.setEDep(delta4[Acts::eEnergy] / Acts::UnitConstants::GeV);
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
