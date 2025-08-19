// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

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
  const std::vector<ScoreMonitor>& scoreMonitor,
  const std::string& monitorFilePath) const {
std::string headers;
headers = "ptScore,";
for (std::size_t i = 0; i < scoreMonitor[0].detectorHitScore.size(); i++) {
  headers += "detectorHitScore" + std::to_string(i) + ",";
}
for (std::size_t i = 0; i < scoreMonitor[0].detectorHoleScore.size(); i++) {
  headers += "detectorHoleScore" + std::to_string(i) + ",";
}
for (std::size_t i = 0; i < scoreMonitor[0].detectorOutlierScore.size();
     i++) {
  headers += "detectorOutlierScore" + std::to_string(i) + ",";
}
for (std::size_t i = 0; i < scoreMonitor[0].detectorOtherScore.size(); i++) {
  headers += "detectorOtherScore" + std::to_string(i) + ",";
}
headers += "chi2Score,";
for (std::size_t i = 0; i < scoreMonitor[0].optionalScore.size(); i++) {
  headers += "optionalScore" + std::to_string(i) + ",";
}
headers += "totalScore\n";

std::ofstream monitorFile;
monitorFile.open(monitorFilePath, std::ios::app);
monitorFile << headers;

for (const auto& monitor : scoreMonitor) {
  monitorFile << monitor.ptScore << ",";
  for (const auto& score : monitor.detectorHitScore) {
    monitorFile << score << ",";
  }
  for (const auto& score : monitor.detectorHoleScore) {
    monitorFile << score << ",";
  }
  for (const auto& score : monitor.detectorOutlierScore) {
    monitorFile << score << ",";
  }
  for (const auto& score : monitor.detectorOtherScore) {
    monitorFile << score << ",";
  }
  monitorFile << monitor.chi2Score << ",";
  for (const auto& score : monitor.optionalScore) {
    monitorFile << score << ",";
  }
  monitorFile << monitor.totalScore << "\n";
}

monitorFile.close();
}