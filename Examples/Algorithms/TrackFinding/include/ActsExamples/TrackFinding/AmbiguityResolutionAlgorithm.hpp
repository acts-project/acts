// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

/// Evicts tracks that seem to be duplicated.
///
/// The implementation works as follows:
///  1) Calculate shared hits per track.
///  2) If the maximum shared hits criteria is met, we are done.
///     This is the configurable amount of shared hits we are ok with
///     in our experiment.
///  3) Else, remove the track with the highest relative shared hits (i.e.
///     shared hits / hits).
///  4) Back to square 1.
class AmbiguityResolutionAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input trajectories collection.
    std::string inputTrajectories;
    /// Input track parameters collection.
    std::string inputTrackParameters;
    /// Input track parameters tips w.r.t outputTrajectories.
    std::string inputTrackParametersTips;
    /// Output track parameters collection.
    std::string outputTrackParameters;
    /// Output track parameters tips w.r.t outputTrajectories.
    std::string outputTrackParametersTips;

    /// Maximum amount of shared hits per track.
    std::uint32_t maximumSharedHits = 1;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  AmbiguityResolutionAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
