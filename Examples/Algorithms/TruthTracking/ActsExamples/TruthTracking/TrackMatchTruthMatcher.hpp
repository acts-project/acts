// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

struct AlgorithmContext;

/// Matches tracks to truth particles and vice versa
class TrackMatchTruthMatcher final : public IAlgorithm {
 public:
  struct Config {
    /// Input (fitted) tracks collection
    std::string inputTracksVT;
    /// Input (fitted) tracks collection
    std::string inputTracksMS;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;

    /// Matching ratio for track to particle matching
    double matchingRatio = 0.5;
  };

  TrackMatchTruthMatcher(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracksVT{this, "InputTracksVT"};
  ReadDataHandle<ConstTrackContainer> m_inputTracksMS{this, "InputTracksMS"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMapVT{
      this, "InputMeasurementParticlesMapVT"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMapMS{
      this, "InputMeasurementParticlesMapMS"};
};

}  // namespace ActsExamples
