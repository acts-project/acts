// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMuonGlobalPatternWriter.hpp"

#include "Acts/Utilities/Enumerate.hpp"

#include <ios>
#include <stdexcept>
#include <unordered_map>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;

namespace ActsExamples {

RootMuonGlobalPatternWriter::RootMuonGlobalPatternWriter(const Config& config,
                                                         Logging::Level level)
    : WriterT(config.inputPatterns, "RootMuonGlobalPatternWriter", level),
      m_cfg{config} {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);

  // open root file and create the tree
  m_file.reset(TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str()));
  if (m_file == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_file->cd();
  m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());

  m_tree->Branch("event_id", &m_eventId);
  m_tree->Branch("pattern_sector", &m_sector);
  m_tree->Branch("pattern_theta", &m_theta);
  m_tree->Branch("pattern_phi", &m_phi);
  m_tree->Branch("pattern_nPrecisionLayers", &m_nPrecisionLayers);
  m_tree->Branch("pattern_nTriggerLayers", &m_nTriggerLayers);
  m_tree->Branch("pattern_nPhiLayers", &m_nPhiLayers);
  m_tree->Branch("pattern_meanNormResidual2", &m_meanNormResidual2);
  m_tree->Branch("pattern_nHits", &m_nHits);

  m_tree->Branch("hit_patternIdx", &m_hitPatternIdx);
  m_tree->Branch("hit_station", &m_hitStation);
  m_tree->Branch("hit_geometryId", &m_hitGeometryId);
  m_tree->Branch("hit_muonId", &m_hitMuonId);
  m_tree->Branch("hit_bucketId", &m_hitBucketId);
  m_tree->Branch("hit_indexInBucket", &m_hitIndexInBucket);
}

RootMuonGlobalPatternWriter::~RootMuonGlobalPatternWriter() = default;

ProcessCode RootMuonGlobalPatternWriter::finalize() {
  m_file->cd();
  m_file->Write();
  m_file.reset();
  ACTS_INFO("Wrote muon global patterns to tree '" << m_cfg.treeName << "' in '"
                                                   << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootMuonGlobalPatternWriter::writeT(
    const AlgorithmContext& ctx, const MuonGlobalPatternContainer& patterns) {
  std::lock_guard lock{m_mutex};
  m_eventId = ctx.eventNumber;

  // The patterns refer to the space points of the event. Their place in the
  // input container identifies a hit across runs of the pattern recognition
  const MuonSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);
  std::unordered_map<const MuonSpacePoint*, std::pair<std::size_t, std::size_t>>
      hitIndices{};
  for (const auto& [bucketIdx, bucket] : enumerate(spacePoints)) {
    for (const auto& [spIdx, spacePoint] : enumerate(bucket)) {
      hitIndices[&spacePoint] = std::make_pair(bucketIdx, spIdx);
    }
  }

  for (const auto& [patternIdx, pattern] : enumerate(patterns)) {
    std::uint16_t nHits{0};
    for (const auto& [station, hits] : enumerate(pattern.hitsPerStation)) {
      for (const MuonSpacePoint* hit : hits) {
        const auto itr = hitIndices.find(hit);
        if (itr == hitIndices.end()) {
          ACTS_WARNING("Pattern hit " << *hit << " is not part of collection "
                                      << m_cfg.inputSpacePoints);
          continue;
        }
        m_hitPatternIdx.push_back(static_cast<std::uint16_t>(patternIdx));
        m_hitStation.push_back(static_cast<std::uint8_t>(station));
        m_hitGeometryId.push_back(hit->geometryId().value());
        m_hitMuonId.push_back(hit->id().toInt());
        m_hitBucketId.push_back(static_cast<std::uint16_t>(itr->second.first));
        m_hitIndexInBucket.push_back(
            static_cast<std::uint16_t>(itr->second.second));
        ++nHits;
      }
    }

    m_sector.push_back(static_cast<std::int16_t>(pattern.sector));
    m_theta.push_back(static_cast<float>(pattern.theta));
    m_phi.push_back(static_cast<float>(pattern.phi));
    m_nPrecisionLayers.push_back(
        static_cast<std::uint8_t>(pattern.nPrecisionLayers));
    m_nTriggerLayers.push_back(
        static_cast<std::uint8_t>(pattern.nTriggerLayers));
    m_nPhiLayers.push_back(static_cast<std::uint8_t>(pattern.nPhiLayers));
    m_meanNormResidual2.push_back(
        static_cast<float>(pattern.meanNormResidual2));
    m_nHits.push_back(nHits);
  }
  ACTS_VERBOSE("Dump " << patterns.size() << " patterns with "
                       << m_hitPatternIdx.size() << " hits.");

  m_tree->Fill();

  m_sector.clear();
  m_theta.clear();
  m_phi.clear();
  m_nPrecisionLayers.clear();
  m_nTriggerLayers.clear();
  m_nPhiLayers.clear();
  m_meanNormResidual2.clear();
  m_nHits.clear();

  m_hitPatternIdx.clear();
  m_hitStation.clear();
  m_hitGeometryId.clear();
  m_hitMuonId.clear();
  m_hitBucketId.clear();
  m_hitIndexInBucket.clear();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
