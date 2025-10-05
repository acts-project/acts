// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMuonSpacePointWriter.hpp"

#include "Acts/Utilities/Enumerate.hpp"

#include "TFile.h"
#include "TTree.h"
namespace ActsExamples {
RootMuonSpacePointWriter::RootMuonSpacePointWriter(const Config& config,
                                                   Acts::Logging::Level level)
    : WriterT(config.inputSpacePoints, "RootMuonSpacePointWriter", level),
      m_cfg{config} {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  m_file.reset(TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str()));
  if (m_file == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_file->cd();
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
}

RootMuonSpacePointWriter::~RootMuonSpacePointWriter() = default;
ProcessCode RootMuonSpacePointWriter::finalize() {
  m_file->cd();
  m_file->Write();
  m_file.reset();
  ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                        << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}
ProcessCode RootMuonSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const MuonSpacePointContainer& hits) {
  m_eventId = ctx.eventNumber;
  for (const auto& [counter, bucket] : Acts::enumerate(hits)) {
    for (const MuonSpacePoint& writeMe : bucket) {
      m_bucketId.push_back(counter);

      m_geometryId.push_back(writeMe.geometryId().value());
      m_muonId.push_back(writeMe.id().toInt());
      m_localPositionX.push_back(writeMe.localPosition().x());
      m_localPositionY.push_back(writeMe.localPosition().y());
      m_localPositionZ.push_back(writeMe.localPosition().z());
      m_sensorDirectionX.push_back(writeMe.sensorDirection().x());
      m_sensorDirectionY.push_back(writeMe.sensorDirection().y());
      m_sensorDirectionZ.push_back(writeMe.sensorDirection().z());
      m_toNextSensorX.push_back(writeMe.toNextSensor().x());
      m_toNextSensorY.push_back(writeMe.toNextSensor().y());
      m_toNextSensorZ.push_back(writeMe.toNextSensor().z());
      m_covLoc0.push_back(writeMe.covariance()[0]);
      m_covLoc1.push_back(writeMe.covariance()[0]);
      m_covLocT.push_back(writeMe.covariance()[0]);
      m_driftR.push_back(writeMe.driftRadius());
      m_time.push_back(writeMe.time());
    }
  }
  m_tree->Fill();

  m_geometryId.clear();
  m_bucketId.clear();
  m_muonId.clear();
  m_localPositionX.clear();
  m_localPositionY.clear();
  m_localPositionZ.clear();
  m_sensorDirectionX.clear();
  m_sensorDirectionY.clear();
  m_sensorDirectionZ.clear();
  m_toNextSensorX.clear();
  m_toNextSensorY.clear();
  m_toNextSensorZ.clear();
  m_covLoc0.clear();
  m_covLoc1.clear();
  m_covLocT.clear();
  m_driftR.clear();
  m_time.clear();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
