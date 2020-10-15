// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"

#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "TrackMlData.hpp"

ActsExamples::CsvSimHitWriter::CsvSimHitWriter(
    const ActsExamples::CsvSimHitWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputSimulatedHits, "CsvSimHitWriter", lvl), m_cfg(cfg) {
    // inputSimHits is already checked by base constructor
    if (m_cfg.outputStem.empty()) {
        throw std::invalid_argument("Missing ouput filename stem");
    }
}

ActsExamples::ProcessCode ActsExamples::CsvSimHitWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::SimHitContainer& simHits) {
  // open per-event file for all simhit components
  std::string pathSimHit = perEventFilepath(
                   m_cfg.outputDir, m_cfg.outputStem + ".csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<SimHitData> writerSimHit(pathSimHit,
                                                     m_cfg.outputPrecision);

  // TrackMlData struct
  SimHitData simhit;
  // Write data from internal impl. to output-side struct
  for (const auto& simHit : simHits) {
    Acts::GeometryIdentifier geoId       = simHit.geometryId();
    // local simhit information in global coord.
    const Acts::Vector4D& globalPos4     = simHit.position4()
    const Acts::Vector4& momentum4Before = simhit.momentum4Before();

    simhit.geometry_id = geoId.value();
    simhit.particle_id = simHit.particleId().value();
    // hit position
    simhit.tx = globalPos4.x() / Acts::UnitConstants::mm;
    simhit.ty = globalPos4.y() / Acts::UnitConstants::mm;
    simhit.tz = globalPos4.z() / Acts::UnitConstants::mm;
    simhit.tt = globalPos4.time() / Acts::UnitConstants::ns;
    // particle four-momentum before interaction
    simhit.tpx = momentum4Before.x() / Acts::UnitConstants::GeV;
    simhit.tpy = momentum4Before.y() / Acts::UnitConstants::GeV;
    simhit.tpz = momentum4Before.z() / Acts::UnitConstants::GeV;
    simhit.te  = momentum4Before.w() / Acts::UnitConstants::GeV;
    // particle four-momentum change due to interaction
    const auto delta4 = simHit.momentum4After() - momentum4Before;
    simhit.deltapx = delta4.x() / Acts::UnitConstants::GeV;
    simhit.deltapy = delta4.y() / Acts::UnitConstants::GeV;
    simhit.deltapz = delta4.z() / Acts::UnitConstants::GeV;
    simhit.deltae  = delta4.w() / Acts::UnitConstants::GeV;
    // TODO write hit index along the particle trajectory
    simhit.index = simHit.index();
    writerSimHit.append(simhit);
    } // end simHit loop

  }

  return ActsExamples::ProcessCode::SUCCESS;
}
