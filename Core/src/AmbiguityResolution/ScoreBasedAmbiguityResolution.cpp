// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include "TFile.h"
#include "TTree.h"

bool Acts::ScoreBasedAmbiguityResolution::etaBasedCuts(
    const DetectorConfig& detector, const TrackFeatures& trackFeatures,
    const double& eta) const {
  const auto& etaBins = detector.etaBins;

  auto it = std::ranges::upper_bound(etaBins, eta);
  if (it == etaBins.begin() || it == etaBins.end()) {
    return false;  // eta out of range
  }
  std::size_t etaBin = std::distance(etaBins.begin(), it) - 1;
  return (trackFeatures.nHits < detector.minHitsPerEta[etaBin] ||
          trackFeatures.nHoles > detector.maxHolesPerEta[etaBin] ||
          trackFeatures.nOutliers > detector.maxOutliersPerEta[etaBin] ||
          trackFeatures.nSharedHits > detector.maxSharedHitsPerEta[etaBin]);
}

void Acts::ScoreBasedAmbiguityResolution::saveScoreMonitor(
    const std::vector<ScoreMonitor>& scoreMonitor) const {
  // Open ROOT file for writing
  TFile* file = TFile::Open(m_monitorFilePath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create ROOT file: " + m_monitorFilePath);
  }

  // Create a TTree
  TTree* tree = new TTree("ScoreMonitorTree", "Score Monitor Data");

  // Variables for branches
  std::size_t ptScore;
  std::size_t chi2Score;
  std::size_t totalScore;

  std::vector<std::size_t> detectorHitScore;
  std::vector<std::size_t> detectorHoleScore;
  std::vector<std::size_t> detectorOutlierScore;
  std::vector<std::size_t> detectorOtherScore;
  std::vector<std::size_t> optionalScore;

  // Create branches
  tree->Branch("ptScore", &ptScore);
  tree->Branch("chi2Score", &chi2Score);
  tree->Branch("totalScore", &totalScore);

  tree->Branch("detectorHitScore", &detectorHitScore);
  tree->Branch("detectorHoleScore", &detectorHoleScore);
  tree->Branch("detectorOutlierScore", &detectorOutlierScore);
  tree->Branch("detectorOtherScore", &detectorOtherScore);
  tree->Branch("optionalScore", &optionalScore);

  // Fill the tree
  for (const ScoreMonitor& monitor : scoreMonitor) {
    ptScore = monitor.ptScore;
    chi2Score = monitor.chi2Score;
    totalScore = monitor.totalScore;

    detectorHitScore = monitor.detectorHitScore;
    detectorHoleScore = monitor.detectorHoleScore;
    detectorOutlierScore = monitor.detectorOutlierScore;
    detectorOtherScore = monitor.detectorOtherScore;
    optionalScore = monitor.optionalScore;

    tree->Fill();
  }

  // Write tree to file
  file->Write();
  file->Close();
  delete file;
}
