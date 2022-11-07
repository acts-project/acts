// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/KDTree.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <optional>
#include <string>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from space points.
class SeedingOrthogonalAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    std::vector<std::string> inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output proto track collection.
    std::string outputProtoTracks;

    Acts::SeedFilterConfig seedFilterConfig;
    Acts::SeedFinderOrthogonalConfig<SimSpacePoint> seedFinderConfig;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingOrthogonalAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const final override;

  /// Const access to the config
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
