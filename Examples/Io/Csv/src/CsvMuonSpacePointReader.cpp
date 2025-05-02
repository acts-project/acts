// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMuonSpacePointReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

#include "CsvOutputData.hpp"

namespace ActsExamples {
CsvMuonSpacePointReader::CsvMuonSpacePointReader(const Config& config,
                                                 Acts::Logging::Level level)
    : m_cfg{config},
      m_eventsRange{
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")},
      m_logger{Acts::getDefaultLogger("CsvMuonSpacePointReader", level)} {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point output collection");
  }

  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
}

std::string CsvMuonSpacePointReader::CsvMuonSpacePointReader::name() const {
  return "CsvMuonSpacePointReader";
}

std::pair<std::size_t, std::size_t> CsvMuonSpacePointReader::availableEvents()
    const {
  return m_eventsRange;
}

ProcessCode CsvMuonSpacePointReader::read(const AlgorithmContext& ctx) {
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);

  NamedTupleCsvReader<MuonSpacePointData> reader(path);

  MuonSpacePointData data{};

  MuonSpacePointContainer spacePoints{};

  using MuonId = MuonSpacePoint::MuonId;
  int lastBucketId{-1};
  while (reader.read(data)) {
    /// Decode the stationName, sector & side of the chamber
    MuonId::StationName stName{
        static_cast<MuonId::StationName>(0x0FF & data.sectorId)};
    MuonId::DetSide side{
        static_cast<MuonId::DetSide>(0x0FF & (data.sectorId >> 8))};
    const auto sector = static_cast<int>(0x0FF & (data.sectorId >> 16));
    MuonId::TechField tech{
        static_cast<MuonId::TechField>(0x0FF & (data.sectorId >> 24))};
    MuonId id{};
    id.setChamber(stName, side, sector, tech);
    id.setLayAndCh(data.gasGap, data.primaryCh);
    id.setCoordFlags(data.measuresEta, data.measuresPhi);
    /// Start a new bucket if the sector id or the bucket Id are different
    if (spacePoints.empty() ||
        !spacePoints.back().back().id().sameStation(id) ||
        lastBucketId != data.bucketId) {
      spacePoints.emplace_back();
      lastBucketId = data.bucketId;
    }
    MuonSpacePoint& newSpacePoint{spacePoints.back().emplace_back()};

    newSpacePoint.setId(id);

    newSpacePoint.defineCoordinates(
        Acts::Vector3{data.locPositionX, data.locPositionY, data.locPositionZ},
        Acts::Vector3{data.locSensorDirX, data.locSensorDirY,
                      data.locSensorDirZ});
    newSpacePoint.defineNormal(Acts::Vector3{
        data.locPlaneNormX, data.locPlaneNormY, data.locPlaneNormZ});
    newSpacePoint.setRadius(data.driftR);

    newSpacePoint.setSpatialCov(data.covXX, data.covXY, data.covYX, data.covYY);
    ACTS_VERBOSE("New spacepoint loaded: " << newSpacePoint);
  }

  // write the ordered data to the EventStore (according to geometry_id).
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
