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
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>
#include <string>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvSpacePointReader::CsvSpacePointReader(const Config& cfg,
                                         Acts::Logging::Level lvl)
    : m_cfg{cfg} {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.inputCollection != "pixel" && m_cfg.inputCollection != "strip" &&
      m_cfg.inputCollection != "overlap") {
    throw std::invalid_argument("Invalid input collection " +
                                m_cfg.inputCollection);
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
  SpacePointContainer spacePoints(
      SpacePointColumns::SourceLinks | SpacePointColumns::X |
      SpacePointColumns::Y | SpacePointColumns::Z |
      SpacePointColumns::VarianceR | SpacePointColumns::VarianceZ |
      SpacePointColumns::Strip);

  const auto& filename = m_cfg.inputCollection.empty()
                             ? m_cfg.inputStem
                             : m_cfg.inputStem + '_' + m_cfg.inputCollection;
  const auto& path =
      perEventFilepath(m_cfg.inputDir, filename + ".csv", ctx.eventNumber);

  NamedTupleCsvReader<SpacePointData> reader(path);
  SpacePointData data;

  while (reader.read(data)) {
    auto sp = spacePoints.createSpacePoint();
    sp.assignSourceLinks(std::array{Acts::SourceLink(data.measurement_id)});
    sp.x() = data.sp_x;
    sp.y() = data.sp_y;
    sp.z() = data.sp_z;
    sp.varianceR() = data.sp_covr;
    sp.varianceZ() = data.sp_covz;

    if (m_cfg.extendCollection) {
      const Acts::Vector3 topStripVector =
          data.sp_topStripDirection * 2 * data.sp_topHalfStripLength;
      const Acts::Vector3 bottomStripVector =
          data.sp_bottomStripDirection * 2 * data.sp_bottomHalfStripLength;
      const Acts::Vector3 stripCenterDistance = data.sp_stripCenterDistance;
      const Acts::Vector3 topStripCenter = data.sp_topStripCenterPosition;

      Eigen::Map<Eigen::Vector3f>(sp.topStripVector().data()) =
          topStripVector.cast<float>();
      Eigen::Map<Eigen::Vector3f>(sp.bottomStripVector().data()) =
          bottomStripVector.cast<float>();
      Eigen::Map<Eigen::Vector3f>(sp.stripCenterDistance().data()) =
          stripCenterDistance.cast<float>();
      Eigen::Map<Eigen::Vector3f>(sp.topStripCenter().data()) =
          topStripCenter.cast<float>();
    }
  }

  ACTS_DEBUG("Created " << spacePoints.size() << " " << m_cfg.inputCollection
                        << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
