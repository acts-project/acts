// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <limits>
#include <string>

namespace ActsExamples {

/// Select tracks by applying some selection cuts.
class TrackModifier final : public IAlgorithm {
 public:
  struct Config {
    /// Optional. Input trajectories container. Mutually exclusive with track
    /// parameters input.
    std::string inputTrajectories;
    /// Optional. Input track parameters collection. Mutually exclusive with
    /// trajectories input.
    std::string inputTrackParameters;
    /// Optional. Output trajectories container. Will only be set if
    /// trajectories input was set
    std::string outputTrajectories;
    /// Optional. Output track parameters collection. Will only be set if track
    /// parameters input was set.
    std::string outputTrackParameters;

    /// When turned on, only keed the diagonal of the cov matrix.
    bool dropCovariance{false};
    /// Scale cov matrix;
    double covScale{1};
    /// Remove all time components
    bool killTime{false};
  };

  TrackModifier(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};
  ReadDataHandle<TrajectoriesContainer> m_inputTrajectories{
      this, "InputTrajectories"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
  WriteDataHandle<TrajectoriesContainer> m_outputTrajectories{
      this, "OutputTrajectories"};
};

}  // namespace ActsExamples
