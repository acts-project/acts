// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <limits>
#include <string>

namespace ActsExamples {
struct AlgorithmContext;

/// Select tracks by applying some selection cuts.
class TrackParameterSelector final : public IAlgorithm {
 public:
  struct Config {
    /// Input track parameters collection
    std::string inputTrackParameters;
    /// Output track parameters collection.
    std::string outputTrackParameters;

    // Minimum/maximum local positions.
    long double loc0Min = -std::numeric_limits<long double>::infinity();
    long double loc0Max = std::numeric_limits<long double>::infinity();
    long double loc1Min = -std::numeric_limits<long double>::infinity();
    long double loc1Max = std::numeric_limits<long double>::infinity();
    // Minimum/maximum track time.
    long double timeMin = -std::numeric_limits<long double>::infinity();
    long double timeMax = std::numeric_limits<long double>::infinity();
    // Direction cuts.
    long double phiMin = -std::numeric_limits<long double>::infinity();
    long double phiMax = std::numeric_limits<long double>::infinity();
    long double etaMin = -std::numeric_limits<long double>::infinity();
    long double etaMax = std::numeric_limits<long double>::infinity();
    long double absEtaMin = 0.0;
    long double absEtaMax = std::numeric_limits<long double>::infinity();
    // Momentum cuts.
    long double ptMin = 0.0;
    long double ptMax = std::numeric_limits<long double>::infinity();
  };

  TrackParameterSelector(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};
  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
};

}  // namespace ActsExamples
