// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootCudaMuonSpacePointReader.hpp"

//#ifdef ACTS_ENABLE_CUDA

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/Io/Root/RootUtility.hpp"

#include <stdexcept>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsExamples {

RootCudaMuonSpacePointReader::RootCudaMuonSpacePointReader(const Config& config, 
                                                           Acts::Logging::Level level)
    : m_cfg{config},
      m_logger{getDefaultLogger("RootCudaMuonSpacePointReader", level)} {
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument(
        "RootCudaMuonSpacePointReader() - Space points must not be empty");
  }

  m_outputContainer.initialize(m_cfg.outputSpacePoints);

  if (!m_file || m_file->IsZombie()) {
    throw std::invalid_argument(
        "RootCudaMuonSpacePointReader() - Failed to open '" + m_cfg.filePath +
        "'");
  }

  if (m_reader.GetTree() == nullptr) {
    throw std::invalid_argument(
        "RootCudaMuonSpacePointReader() - Failed to load TTree '" +
        m_cfg.treeName + "'");
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

RootCudaMuonSpacePointReader::~RootCudaMuonSpacePointReader() = default;

std::pair<std::size_t, std::size_t>
RootCudaMuonSpacePointReader::availableEvents() const {
  return std::make_pair(0u, m_reader.GetEntries());
}

ProcessCode RootCudaMuonSpacePointReader::read(
    const AlgorithmContext& context) {
  std::lock_guard guard{m_mutex};

  const auto entry = m_eventRanges.at(context.eventNumber);
  m_reader.SetEntry(entry);

  CudaMuonSpacePointContainer outSpacePoints{m_bucketId->size()};

  std::size_t bucketStart = 0;
  std::uint16_t currentBucket = 0;
  bool hasOpenBucket = false;

  for (std::size_t spIdx = 0; spIdx < m_bucketId->size(); ++spIdx) {
    const auto bucketIdx = m_bucketId->at(spIdx);

    if (!hasOpenBucket) {
      bucketStart = spIdx;
      currentBucket = bucketIdx;
      hasOpenBucket = true;
    } else if (bucketIdx != currentBucket) {
      outSpacePoints.addBucket(bucketStart, spIdx);
      bucketStart = spIdx;
      currentBucket = bucketIdx;
    }

    outSpacePoints.setGeometryId(spIdx, m_geometryId->at(spIdx));
    outSpacePoints.setId(spIdx, m_muonId->at(spIdx));

    Vector3 position{m_localPositionX->at(spIdx), m_localPositionY->at(spIdx),
                     m_localPositionZ->at(spIdx)};

    Vector3 sensorDir{makeDirectionFromPhiTheta<double>(
        m_sensorDirectionPhi->at(spIdx) * 1._degree,
        m_sensorDirectionTheta->at(spIdx) * 1._degree)};

    Vector3 toNext{makeDirectionFromPhiTheta<double>(
        m_toNextSensorPhi->at(spIdx) * 1._degree,
        m_toNextSensorTheta->at(spIdx) * 1._degree)};

    outSpacePoints.defineCoordinates(spIdx, position, sensorDir, toNext);

    outSpacePoints.setRadius(spIdx, m_driftR->at(spIdx));
    outSpacePoints.setTime(spIdx, m_time->at(spIdx));
    outSpacePoints.setCovariance(spIdx, m_covLoc0->at(spIdx),
                                 m_covLoc1->at(spIdx), m_covT->at(spIdx));

    ACTS_VERBOSE("Loaded CUDA muon space point " << spIdx);
  }

  if (hasOpenBucket) {
    outSpacePoints.addBucket(bucketStart, m_bucketId->size());
  }

  m_outputContainer(context, std::move(outSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

#endif
