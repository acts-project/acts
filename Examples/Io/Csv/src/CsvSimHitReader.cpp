// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvSimHitReader::CsvSimHitReader(
    const ActsExamples::CsvSimHitReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvSimHitReader", level)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }

  m_outputSimHits.initialize(m_cfg.outputSimHits);
}

std::string ActsExamples::CsvSimHitReader::CsvSimHitReader::name() const {
  return "CsvSimHitReader";
}

std::pair<size_t, size_t> ActsExamples::CsvSimHitReader::availableEvents()
    const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::CsvSimHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);

  dfe::NamedTupleCsvReader<SimHitData> reader(path);

  SimHitContainer::sequence_type unordered;
  SimHitData data;

  while (reader.read(data)) {
    const auto geometryId = Acts::GeometryIdentifier(data.geometry_id);
    // TODO validate geo id consistency
    const auto particleId = ActsFatras::Barcode(data.particle_id);

    ActsFatras::Hit::Vector4 pos4{
        data.tx * Acts::UnitConstants::mm,
        data.ty * Acts::UnitConstants::mm,
        data.tz * Acts::UnitConstants::mm,
        data.tt * Acts::UnitConstants::ns,
    };
    ActsFatras::Hit::Vector4 mom4{
        data.tpx * Acts::UnitConstants::GeV,
        data.tpy * Acts::UnitConstants::GeV,
        data.tpz * Acts::UnitConstants::GeV,
        data.te * Acts::UnitConstants::GeV,
    };
    ActsFatras::Hit::Vector4 delta4{
        data.deltapx * Acts::UnitConstants::GeV,
        data.deltapy * Acts::UnitConstants::GeV,
        data.deltapz * Acts::UnitConstants::GeV,
        data.deltae * Acts::UnitConstants::GeV,
    };

    ActsFatras::Hit hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                        data.index);
    unordered.push_back(std::move(hit));
  }

  // write the ordered data to the EventStore (according to geometry_id).
  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());
  m_outputSimHits(ctx, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
