// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootCudaMuonSpacePointWriter.hpp"

#ifdef ACTS_ENABLE_CUDA

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <ios>
#include <limits>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace Acts::VectorHelpers;

namespace {

/// Pushes value back to a vector and applies a static cast.
/// @param vec Vector into which the value is pushed back.
/// @param val Value to push.
template <typename T, typename T1>
void castPush(std::vector<T>& vec, const T1& val) {
  const T castedVal = static_cast<T>(val);
  vec.push_back(castedVal);
}

/// Converts an angle in radians into an angle in degree.
constexpr double inDeg(const double radians) {
  if (Acts::abs(radians) < std::numeric_limits<float>::epsilon()) {
    return 0.;
  }

  using namespace Acts::UnitLiterals;
  return radians / 1._degree;
}

}  // namespace

namespace ActsExamples {

RootCudaMuonSpacePointWriter::RootCudaMuonSpacePointWriter(
    const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputSpacePoints, "RootCudaMuonSpacePointWriter", level),
      m_cfg{config} {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument(
        "RootCudaMuonSpacePointWriter - Missing file path");
  }

  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument(
        "RootCudaMuonSpacePointWriter - Missing tree name");
  }

  m_file.reset(TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str()));
  if (m_file == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_file->cd();
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  m_tree->Branch("event_id", &m_eventId);

  m_tree->Branch("spacePoint_bucketId", &m_bucketId);
  m_tree->Branch("spacePoint_geometryId", &m_geometryId);
  m_tree->Branch("spacePoint_muonId", &m_muonId);

  m_tree->Branch("spacePoint_localPosX", &m_localPositionX);
  m_tree->Branch("spacePoint_localPosY", &m_localPositionY);
  m_tree->Branch("spacePoint_localPosZ", &m_localPositionZ);

  m_tree->Branch("spacePoint_sensorDirTheta", &m_sensorDirectionTheta);
  m_tree->Branch("spacePoint_sensorDirPhi", &m_sensorDirectionPhi);

  m_tree->Branch("spacePoint_toNextDirTheta", &m_toNextSensorTheta);
  m_tree->Branch("spacePoint_toNextDirPhi", &m_toNextSensorPhi);

  m_tree->Branch("spacePoint_covLoc0", &m_covLoc0);
  m_tree->Branch("spacePoint_covLoc1", &m_covLoc1);
  m_tree->Branch("spacePoint_covT", &m_covT);

  m_tree->Branch("spacePoint_driftRadius", &m_driftR);
  m_tree->Branch("spacePoint_time", &m_time);
}

RootCudaMuonSpacePointWriter::~RootCudaMuonSpacePointWriter() = default;

ProcessCode RootCudaMuonSpacePointWriter::finalize() {
  m_file->cd();
  m_file->Write();
  m_file.reset();

  ACTS_INFO("Wrote CUDA muon space points to tree '"
            << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootCudaMuonSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const CudaMuonSpacePointContainer& hits) {
  std::lock_guard lock{m_mutex};

  m_eventId = ctx.eventNumber;

  for (std::size_t bucketIdx = 0; bucketIdx < hits.bucketCount(); ++bucketIdx) {
    for (std::size_t spIdx = hits.bucketStart(bucketIdx);
         spIdx < hits.bucketEnd(bucketIdx); ++spIdx) {
      auto writeMe = hits[spIdx];

      castPush(m_bucketId, bucketIdx);
      castPush(m_geometryId, writeMe->geometryId().value());
      castPush(m_muonId, writeMe->id().toInt());

      castPush(m_localPositionX, writeMe->localPosition().x());
      castPush(m_localPositionY, writeMe->localPosition().y());
      castPush(m_localPositionZ, writeMe->localPosition().z());

      castPush(m_sensorDirectionTheta,
               inDeg(theta(writeMe->sensorDirection())));
      castPush(m_sensorDirectionPhi, inDeg(phi(writeMe->sensorDirection())));

      castPush(m_toNextSensorTheta, inDeg(theta(writeMe->toNextSensor())));
      castPush(m_toNextSensorPhi, inDeg(phi(writeMe->toNextSensor())));

      const auto& cov = writeMe->covariance();
      castPush(m_covLoc0, cov[0]);
      castPush(m_covLoc1, cov[1]);
      castPush(m_covT, cov[2]);

      castPush(m_driftR, writeMe->driftRadius());
      castPush(m_time, writeMe->time());

      ACTS_VERBOSE("Dump CUDA muon space point " << spIdx);
    }
  }

  m_tree->Fill();

  m_geometryId.clear();
  m_bucketId.clear();
  m_muonId.clear();

  m_localPositionX.clear();
  m_localPositionY.clear();
  m_localPositionZ.clear();

  m_sensorDirectionTheta.clear();
  m_sensorDirectionPhi.clear();

  m_toNextSensorTheta.clear();
  m_toNextSensorPhi.clear();

  m_covLoc0.clear();
  m_covLoc1.clear();
  m_covT.clear();

  m_driftR.clear();
  m_time.clear();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

#endif
