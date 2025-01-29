// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

std::size_t Acts::ScoreBasedAmbiguityResolution::getValueAtEta(
    std::vector<std::size_t> cuts, std::size_t etaBinSize,
    std::size_t binIndex) {
  if (cuts.size() == etaBinSize - 1) {
    return cuts[binIndex];
  }

  else if (cuts.size() == 1) {
    return cuts[0];
  }

  else {
    throw std::invalid_argument("Invalid cuts size. Expected 1 or " +
                                std::to_string(etaBinSize - 1) + ", got " +
                                std::to_string(cuts.size()));
  }
}

bool Acts::ScoreBasedAmbiguityResolution::etaBasedCuts(
    const DetectorConfig& detector, const TrackFeatures& trackFeatures,
    const double& eta) const {
  std::vector<double> etaBins = detector.etaBins;

  auto it = std::upper_bound(etaBins.begin(), etaBins.end(), eta);
  if (it == etaBins.begin() || it == etaBins.end()) {
    return false;  // eta out of range
  }
  std::size_t etaBin = std::distance(etaBins.begin(), it) - 1;
  return (trackFeatures.nHits <
          getValueAtEta(detector.minHitsPerEta, etaBins.size(), etaBin)) ||
         (trackFeatures.nHoles >
          getValueAtEta(detector.maxHolesPerEta, etaBins.size(), etaBin)) ||
         (trackFeatures.nOutliers >
          getValueAtEta(detector.maxOutliersPerEta, etaBins.size(), etaBin));
}
