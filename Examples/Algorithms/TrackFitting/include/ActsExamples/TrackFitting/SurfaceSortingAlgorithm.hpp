// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <map>
#include <memory>
#include <vector>

namespace ActsExamples {

using TrackHitList = std::map<const double, const Index>;

class SurfaceSortingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input proto track collection
    std::string inputProtoTracks;
    /// Input simulated hit collection
    std::string inputSimHits;
    /// Input measurement to simulated hit map for truth position
    std::string inputMeasurementSimHitsMap;
    /// Output proto track collection
    std::string outputProtoTracks;
  };

  SurfaceSortingAlgorithm(Config cfg, Acts::Logging::Level level);

  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
