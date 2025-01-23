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
    std::cerr
        << "The size of minHits is not equal to 1 or the number of eta bins"
        << std::endl;
    return 0;
  }
}

bool Acts::ScoreBasedAmbiguityResolution::etaBasedCuts(
    const DetectorConfig& detector, const TrackFeatures& trackFeatures,
    const double& eta) const {
  std::vector<double> etaBins = detector.etaBins;
  int etaBin = 0;

  bool cutApplied = false;

  for (std::size_t i = 0; i < etaBins.size() - 1; ++i) {
    if (eta >= etaBins[i] && eta < etaBins[i + 1]) {
      etaBin = i;
      break;
    }
  }

  cutApplied = (trackFeatures.nHits <= getValueAtEta(detector.minHitsPerEta,
                                                     etaBins.size(), etaBin)) ||
               cutApplied;
  cutApplied =
      (trackFeatures.nHoles >=
       getValueAtEta(detector.maxHolesPerEta, etaBins.size(), etaBin)) ||
      cutApplied;
  cutApplied =
      (trackFeatures.nOutliers >=
       getValueAtEta(detector.maxOutliersPerEta, etaBins.size(), etaBin)) ||
      cutApplied;

  return cutApplied;
}
