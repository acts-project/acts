// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <optional>
#include <stdexcept>
#include <string>

#include <boost/container/static_vector.hpp>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvSpacePointReader::CsvSpacePointReader(const Config& cfg,
                                         Acts::Logging::Level lvl) {
  m_cfg = cfg;
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  auto& filename = m_cfg.inputCollection.empty()
                       ? cfg.inputStem
                       : cfg.inputStem + '_' + cfg.inputCollection;
  m_eventsRange = determineEventFilesRange(cfg.inputDir, filename + ".csv");
  m_logger = Acts::getDefaultLogger("CsvSpacePointReader", lvl);

  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
}

std::string CsvSpacePointReader::CsvSpacePointReader::name() const {
  return "CsvSpacePointReader";
}

std::pair<std::size_t, std::size_t> CsvSpacePointReader::availableEvents()
    const {
  return m_eventsRange;
}

ProcessCode CsvSpacePointReader::read(const AlgorithmContext& ctx) {
  SimSpacePointContainer spacePoints;

  const auto& filename = m_cfg.inputCollection.empty()
                             ? m_cfg.inputStem
                             : m_cfg.inputStem + '_' + m_cfg.inputCollection;
  const auto& path =
      perEventFilepath(m_cfg.inputDir, filename + ".csv", ctx.eventNumber);

  BoostDescribeCsvReader<SpacePointData> reader(path);
  SpacePointData data;

  while (reader.read(data)) {
    Acts::Vector3 globalPos(data.sp_x, data.sp_y, data.sp_z);

    if (m_cfg.inputCollection == "pixel" || m_cfg.inputCollection == "strip" ||
        m_cfg.inputCollection == "overlap") {
      boost::container::static_vector<Acts::SourceLink, 2> sLinks;
      // auto sp = SimSpacePoint(globalPos, data.sp_covr, data.sp_covz, sLinks);

      if (m_cfg.extendCollection) {
        Acts::Vector3 topStripDirection(data.sp_topStripDirection_0,
                                        data.sp_topStripDirection_1,
                                        data.sp_topStripDirection_2);
        Acts::Vector3 bottomStripDirection(data.sp_bottomStripDirection_0,
                                           data.sp_bottomStripDirection_1,
                                           data.sp_bottomStripDirection_2);
        Acts::Vector3 stripCenterDistance(data.sp_stripCenterDistance_0,
                                          data.sp_stripCenterDistance_1,
                                          data.sp_stripCenterDistance_2);
        Acts::Vector3 topStripCenterPosition(data.sp_topStripCenterPosition_0,
                                             data.sp_topStripCenterPosition_1,
                                             data.sp_topStripCenterPosition_2);

        // TODO time
        spacePoints.emplace_back(
            globalPos, std::nullopt, data.sp_covr, data.sp_covz, std::nullopt,
            sLinks, data.sp_topHalfStripLength, data.sp_bottomHalfStripLength,
            topStripDirection, bottomStripDirection, stripCenterDistance,
            topStripCenterPosition);
      } else {
        // TODO time
        spacePoints.emplace_back(globalPos, std::nullopt, data.sp_covr,
                                 data.sp_covz, std::nullopt, sLinks);
      }
    } else {
      ACTS_ERROR("Invalid space point type " << m_cfg.inputStem);
      return ProcessCode::ABORT;
    }
  }

  ACTS_DEBUG("Created " << spacePoints.size() << " " << m_cfg.inputCollection
                        << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
