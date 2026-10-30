// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFinding/SeedParamsEstimator.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace ActsExamples {

/// Estimate track parameters from the space points of a seed.
///
/// Seeds without an estimate are dropped. Optional input seeds and proto
/// tracks are written back out filtered the same way.
class TrackParamsEstimationAlgorithm final : public IAlgorithm {
 public:
  struct Config : SeedParamsEstimator::Config {
    /// Input seeds collection.
    std::string inputSeeds;
    /// Input proto tracks (optional).
    std::optional<std::string> inputProtoTracks;
    /// Input particle hypothesis (optional). If not given, the static particle
    /// hypothesis from the config is used.
    std::optional<std::string> inputParticleHypotheses;
    /// Output estimated track parameters collection.
    std::string outputTrackParameters;
    /// Output seed collection - only seeds with successful parameter estimation
    /// are propagated (optional)
    std::optional<std::string> outputSeeds;
    /// Output proto track collection - only tracks with successful parameter
    /// estimation are propagated (optional)
    std::optional<std::string> outputProtoTracks;

    /// Input free track parameters per seed, used instead of the fit to the
    /// space points (optional).
    std::optional<std::string> inputFreeParameters;
    /// Give the parameters on the surface of the last space point of the seed
    /// instead of the first, for a track finder that starts there.
    bool parametersOnLastSpacePoint = false;
  };

  /// Construct the track parameters making algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit TrackParamsEstimationAlgorithm(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the track parameters making algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  SeedParamsEstimator m_estimator;

  ReadDataHandle<SeedContainer> m_inputSeeds{this, "InputSeeds"};
  ReadDataHandle<ProtoTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<std::vector<Acts::ParticleHypothesis>>
      m_inputParticleHypotheses{this, "InputParticleHypotheses"};

  ReadDataHandle<std::vector<Acts::FreeVector>> m_inputFreeParameters{
      this, "InputFreeParameters"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
  WriteDataHandle<SeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ProtoTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
