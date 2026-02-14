// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMuonSpacePointReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"
using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

RootMuonSpacePointReader::RootMuonSpacePointReader(const Config& config,
                                                   Acts::Logging::Level level)
    : m_cfg{config},
      m_logger{getDefaultLogger("RootMuonSpacePointReader", level)} {
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument(
        "RootMuonSpacePointReader() - Space points must not be empty");
  }
  m_outputContainer.initialize(m_cfg.outputSpacePoints);
  if (!m_file || m_file->IsZombie()) {
    throw std::invalid_argument(
        "RootMuonSpacePointReader() - Failed to open '" + m_cfg.filePath + "'");
  }
  if (m_reader.GetTree() == nullptr) {
    throw std::invalid_argument("Failed to load TTree '" + m_cfg.treeName +
                                "'");
  }

  // Sort the entry numbers of the events
  {
    // necessary to guarantee that m_inputChain->GetV1() is valid for the
    // entire range
    m_eventRanges.resize(m_reader.GetEntries());
    m_reader.GetTree()->SetEstimate(m_eventRanges.size() + 1);

    m_reader.GetTree()->Draw("event_id", "", "goff");
    const auto nEntries = static_cast<std::uint32_t>(m_eventRanges.size());
    RootUtility::stableSort(nEntries, m_reader.GetTree()->GetV1(),
                            m_eventRanges.data(), false);
  }
}

RootMuonSpacePointReader::~RootMuonSpacePointReader() = default;

std::pair<std::size_t, std::size_t> RootMuonSpacePointReader::availableEvents()
    const {
  return std::make_pair(0u, m_reader.GetEntries());
}

ProcessCode RootMuonSpacePointReader::read(const AlgorithmContext& context) {
  std::lock_guard guard{m_mutex};

  MuonSpacePointContainer outSpacePoints{};

  auto entry = m_eventRanges.at(context.eventNumber);
  m_reader.SetEntry(entry);
  for (std::size_t spIdx = 0; spIdx < m_bucketId->size(); ++spIdx) {
    auto bucketIdx = static_cast<std::size_t>(m_bucketId->at(spIdx));
    // The space point buckets are ordered sequentially
    if (bucketIdx + 1 != outSpacePoints.size()) {
      outSpacePoints.emplace_back();
    }
    MuonSpacePoint& newSp{outSpacePoints.back().emplace_back()};
    newSp.setGeometryId(GeometryIdentifier{m_geometryId->at(spIdx)});
    newSp.setId(MuonSpacePoint::MuonId{m_muonId->at(spIdx)});

    Vector3 position{m_localPositionX->at(spIdx), m_localPositionY->at(spIdx),
                     m_localPositionZ->at(spIdx)};
    Vector3 sensorDir{makeDirectionFromPhiTheta<double>(
        m_sensorDirectionPhi->at(spIdx) * 1._degree,
        m_sensorDirectionTheta->at(spIdx) * 1._degree)};
    Vector3 toNext{makeDirectionFromPhiTheta<double>(
        m_toNextSensorPhi->at(spIdx) * 1._degree,
        m_toNextSensorTheta->at(spIdx) * 1._degree)};

    newSp.defineCoordinates(std::move(position), std::move(sensorDir),
                            std::move(toNext));

    newSp.setRadius(m_driftR->at(spIdx));
    newSp.setTime(m_time->at(spIdx));
    newSp.setCovariance(m_covLoc0->at(spIdx), m_covLoc1->at(spIdx),
                        m_covT->at(spIdx));
    ACTS_VERBOSE("Loaded new space point " << newSp);
  }

  m_outputContainer(context, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
