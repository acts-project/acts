// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"

#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <dfe/dfe_io_dsv.hpp>

#include "TrackMlData.hpp"

ActsExamples::CsvSimHitReader::CsvSimHitReader(
    const ActsExamples::CsvSimHitReader::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg)
      // TODO check that all files (hits,cells,truth) exists
      ,
      m_eventsRange(determineEventFilesRange(cfg.inputDir, "hits.csv")),
      m_logger(Acts::getDefaultLogger("CsvSimHitReader", lvl)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }
}

std::string ActsExamples::CsvSimHitReader::CsvSimHitReader::name()
    const {
  return "CsvSimHitReader";
}

std::pair<size_t, size_t>
ActsExamples::CsvSimHitReader::availableEvents() const {
  return m_eventsRange;
}


ActsExamples::ProcessCode ActsExamples::CsvSimHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  // geometry_id in the files is not required to be neither continuous nor
  // monotonic. internally, we want continous indices within [0,#geoObj)
  // to simplify data handling. to be able to perform this mapping we first
  // read all data into memory before converting to the internal event data
  // types.
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);
  // define all optional columns
  std::vector<std::string> optionalColumns = {
      "geometry_id", "tt",      "te",     "deltapx",
      "deltapy",     "deltapz", "deltae", "index",
  };
  dfe::NamedTupleCsvReader<SimHitData> reader(path, optionalColumns);

  std::vector<ActsFatras::Hit> unordered;
  SimHitData data;

  while (reader.read(data)) {
    const auto geometryId  = Acts::GeometryIdentifier(simhit.geometry_id);
    // TODO validate geo id consistency
    const auto particleId  = ActsFatras::Barcode(simhit.particle_id);

    ActsFatras::Hit::Vector4 pos4{
      simhit.tx * Acts::UnitConstants::mm,
      simhit.ty * Acts::UnitConstants::mm,
      simhit.tz * Acts::UnitConstants::mm,
      simhit.tt * Acts::UnitConstants::ns,
    };
    ActsFatras::Hit::Vector4 mom4{
      simhit.tpx * Acts::UnitConstants::GeV,
      simhit.tpy * Acts::UnitConstants::GeV,
      simhit.tpz * Acts::UnitConstants::GeV,
      simhit.te  * Acts::UnitConstants::GeV,
    };
    ActsFatras::Hit::Vector4 delta4{
      simhit.deltapx * Acts::UnitConstants::GeV,
      simhit.deltapy * Acts::UnitConstants::GeV,
      simhit.deltapz * Acts::UnitConstants::GeV,
      simhit.deltae  * Acts::UnitConstants::GeV,
    };

    ActsFatras::Hit hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                        simhit.index());
    unordered.push_back(std::move(hit));
  }

  // write the ordered data to the EventStore
  SimHitContainer simHits;
  simHits.adopt_sequence(std::move(unordered));
  ctx.eventStore.add(m_cfg.outputSimulatedHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
