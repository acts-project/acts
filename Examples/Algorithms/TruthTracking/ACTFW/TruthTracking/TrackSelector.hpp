// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <string>

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Select tracks by applying some selection cuts.
class TrackSelector final : public BareAlgorithm {
 public:
  struct Config {
    /// The input collection
    std::string input;
    /// The output collection
    std::string output;
    /// Maximum distance from the origin in the transverse plane
    double rhoMax = std::numeric_limits<double>::max();
    /// Maximum absolute distance from the origin along z
    double absZMax = std::numeric_limits<double>::max();
    // Track cuts
    double phiMin = std::numeric_limits<double>::lowest();
    double phiMax = std::numeric_limits<double>::max();
    double etaMin = std::numeric_limits<double>::lowest();
    double etaMax = std::numeric_limits<double>::max();
    double absEtaMin = std::numeric_limits<double>::lowest();
    double absEtaMax = std::numeric_limits<double>::max();
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::max();
    /// Keep neutral particles
    bool keepNeutral = true;
  };

  TrackSelector(const Config& cfg,
                Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
