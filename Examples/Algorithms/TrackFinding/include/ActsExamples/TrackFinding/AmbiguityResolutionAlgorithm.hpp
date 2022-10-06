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

///
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
    /// Output track indices.
    std::string outputTrackIndices;

    ///
    std::uint32_t maximumSharedHits = 1;
  };

  /// 
  AmbiguityResolutionAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// 
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// 
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
