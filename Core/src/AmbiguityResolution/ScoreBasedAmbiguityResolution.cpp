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
  std::vector<double> etaBins = detector.etaBins;

  auto it = std::upper_bound(etaBins.begin(), etaBins.end(), eta);
  if (it == etaBins.begin() || it == etaBins.end()) {
    return false;  // eta out of range
  }
  std::size_t etaBin = std::distance(etaBins.begin(), it) - 1;
  return (trackFeatures.nHits < detector.minHitsPerEta[etaBin] ||
          trackFeatures.nHoles > detector.maxHolesPerEta[etaBin] ||
          trackFeatures.nOutliers > detector.maxOutliersPerEta[etaBin]);
}
