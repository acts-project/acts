// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMuonSegmentReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/MuonSegment.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <bitset>
#include <stdexcept>

#include "CsvOutputData.hpp"

namespace ActsExamples {
CsvMuonSegmentReader::CsvMuonSegmentReader(
    const CsvMuonSegmentReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvMuonSegmentReader", level)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputSegments.empty()) {
    throw std::invalid_argument("Missing segment collection");
  }

  m_outputSegments.initialize(m_cfg.outputSegments);
}

std::string CsvMuonSegmentReader::CsvMuonSegmentReader::name() const {
  return "CsvMuonSegmentReader";
}

std::pair<std::size_t, std::size_t> CsvMuonSegmentReader::availableEvents()
    const {
  return m_eventsRange;
}

ProcessCode CsvMuonSegmentReader::read(const AlgorithmContext& ctx) {
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);

  NamedTupleCsvReader<MuonSegmentData> reader(path);

  MuonSegmentData data{};

  MuonSegmentContainer segments{};

  using MuonId = MuonSegment::MuonId;
  while (reader.read(data)) {
    MuonSegment& newSeg = segments.emplace_back();
    /// Decode the stationName, sector & side of the chamber
    MuonId::StationName stName{
        static_cast<MuonId::StationName>(0x0FF & data.sectorId)};
    MuonId::DetSide side{
        static_cast<MuonId::DetSide>(0x0FF & (data.sectorId >> 8))};
    const auto sector = static_cast<int>(0x0FF & (data.sectorId >> 16));
    MuonId id{};
    id.setChamber(stName, side, sector, MuonId::TechField::UnDef);
    ACTS_VERBOSE("Read in new segment in "
                 << id << ", sector id: " << data.sectorId << ", "
                 << std::bitset<32>(data.sectorId));
    newSeg.setId(id);
    newSeg.setGlobalCoords(
        Acts::Vector3{data.globalPositionX, data.globalPositionY,
                      data.globalPositionZ},
        Acts::Vector3{data.globalDirectionX, data.globalDirectionY,
                      data.globalDirectionZ});
    newSeg.setLocalCoords(
        Acts::Vector3{data.localPositionX, data.localPositionY,
                      data.localPositionZ},
        Acts::Vector3{data.localDirectionX, data.localDirectionY,
                      data.localDirectionZ});
    newSeg.setTime(data.time, data.timeError);
    newSeg.setFitQuality(data.chiSquared, data.nDoF);
    newSeg.setHitSummary(data.precisionHits, data.trigEtaLayers,
                         data.phiLayers);
  }
  // write the ordered data to the EventStore (according to geometry_id).
  m_outputSegments(ctx, std::move(segments));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
