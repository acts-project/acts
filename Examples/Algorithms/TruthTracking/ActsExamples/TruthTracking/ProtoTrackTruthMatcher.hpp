// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

struct AlgorithmContext;

/// Matches proto track to truth particles and vice versa
class ProtoTrackTruthMatcher final : public IAlgorithm {
 public:
  struct Config {
    /// Input proto tracks collection
    std::string inputProtoTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Output proto track-particle matching.
    std::string outputProtoTrackParticleMatching;
    /// Output particle-proto track matching.
    std::string outputParticleProtoTrackMatching;

    /// Matching ratio for track to particle matching
    double matchingRatio = 0.5;
    /// Whether to use double matching (track to particle and particle to track)
    bool doubleMatching = false;
  };

  ProtoTrackTruthMatcher(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  WriteDataHandle<TrackParticleMatching> m_outputProtoTrackParticleMatching{
      this, "OutputProtoTrackParticleMatching"};
  WriteDataHandle<ParticleTrackMatching> m_outputParticleProtoTrackMatching{
      this, "OutputParticleProtoTrackMatching"};
};

}  // namespace ActsExamples
