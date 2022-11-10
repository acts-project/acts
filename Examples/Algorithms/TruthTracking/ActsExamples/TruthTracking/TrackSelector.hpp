// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <limits>
#include <string>

namespace ActsExamples {

/// Select tracks by applying some selection cuts.
class TrackSelector final : public BareAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection. Mutually exclusive with
    /// trajectories input.
    std::string inputTrackParameters;
    /// Optional. Input trajectories container. Mutually exclusive with track
    /// parameters input.
    std::string inputTrajectories;
    /// Optional. Output track parameters collection. Will only be set if track
    /// parameters input was set.
    std::string outputTrackParameters;
    /// Optional. Output trajectories container. Will only be set if
    /// trajectories input was set
    std::string outputTrajectories;

    // Minimum/maximum local positions.
    double loc0Min = -std::numeric_limits<double>::infinity();
    double loc0Max = std::numeric_limits<double>::infinity();
    double loc1Min = -std::numeric_limits<double>::infinity();
    double loc1Max = std::numeric_limits<double>::infinity();
    // Minimum/maximum track time.
    double timeMin = -std::numeric_limits<double>::infinity();
    double timeMax = std::numeric_limits<double>::infinity();
    // Direction cuts.
    double phiMin = -std::numeric_limits<double>::infinity();
    double phiMax = std::numeric_limits<double>::infinity();
    double etaMin = -std::numeric_limits<double>::infinity();
    double etaMax = std::numeric_limits<double>::infinity();
    double absEtaMin = 0.0;
    double absEtaMax = std::numeric_limits<double>::infinity();
    // Momentum cuts.
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::infinity();
    /// Remove charged particles.
    bool removeCharged = false;
    /// Remove neutral particles.
    bool removeNeutral = false;
  };

  TrackSelector(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
