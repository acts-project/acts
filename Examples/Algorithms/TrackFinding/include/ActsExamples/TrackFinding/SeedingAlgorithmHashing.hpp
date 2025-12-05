// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/detail/CylindricalSpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/SpacePointContainer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Hashing/HashingAlgorithm.hpp"
#include "ActsPlugins/Hashing/HashingAlgorithmConfig.hpp"
#include "ActsPlugins/Hashing/HashingTraining.hpp"
#include "ActsPlugins/Hashing/HashingTrainingConfig.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from space points.
class SeedingAlgorithmHashing final : public IAlgorithm {
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
    /// Output space point buckets.
    std::string outputBuckets;

    Acts::SeedFilterConfig seedFilterConfig;
    Acts::SeedFinderConfig<typename Acts::SpacePointContainer<
        ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
        Acts::detail::RefHolder>::SpacePointProxyType>
        seedFinderConfig;

    Acts::CylindricalSpacePointGridConfig gridConfig;
    Acts::CylindricalSpacePointGridOptions gridOptions;
    Acts::SeedFinderOptions seedFinderOptions;
    ActsPlugins::HashingAlgorithmConfig hashingConfig;
    ActsPlugins::HashingTrainingConfig hashingTrainingConfig;

    // allow for different values of rMax in gridConfig and seedFinderConfig
    bool allowSeparateRMax = false;

    // vector containing the map of z bins in the top and bottom layers
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;
    // number of phiBin neighbors at each side of the current bin that will be
    // used to search for SPs
    int numPhiNeighbors = 1;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithmHashing(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  using SpacePointProxy_t = typename Acts::SpacePointContainer<
      ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
      Acts::detail::RefHolder>::SpacePointProxyType;

  Acts::SeedFinder<SpacePointProxy_t,
                   Acts::CylindricalSpacePointGrid<SpacePointProxy_t>>
      m_seedFinder;
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_bottomBinFinder{nullptr};
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_topBinFinder{nullptr};

  Config m_cfg;

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<std::vector<SimSpacePointContainer>> m_outputBuckets{
      this, "OutputBuckets"};
  ActsPlugins::HashingAlgorithm<const SimSpacePoint*,
                                std::vector<const SimSpacePoint*>>
      m_hashing;
  ActsPlugins::HashingTrainingAlgorithm<std::vector<const SimSpacePoint*>>
      m_hashingTraining;
};

}  // namespace ActsExamples
