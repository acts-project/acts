// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData2/SeedContainer2.hpp"
#include "Acts/EventData2/SpacePointContainer2.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding2/IExperimentCuts2.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace Acts {

class SeedFilter2 final {
 public:
  struct DerivedConfig;

  /// @brief Structure that holds configuration parameters for the seed filter algorithm
  struct Config {
    /// Allowed difference in curvature (inverted seed radii) between two
    /// compatible seeds
    float deltaInvHelixDiameter = 0.00003 * 1. / UnitConstants::mm;
    /// Minimum distance between compatible outer space-points to be considered.
    /// This is used to avoid counting space-points from the same layer
    float deltaRMin = 5. * UnitConstants::mm;
    /// Seed weight/score is increased by this value if a compatible seed has
    /// been found. This is the c1 factor in the seed score calculation (w = c1
    /// * Nt - c2 * d0 - c3 * z0)
    float compatSeedWeight = 200.;
    /// The transverse impact parameters (d0) is multiplied by this factor and
    /// subtracted from weight. This is the c2 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float impactWeightFactor = 1.;
    /// The logitudinal impact parameters (z0) is multiplied by this factor and
    /// subtracted from weight. This is the c3 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float zOriginWeightFactor = 1.;
    /// Maximum number (minus one) of accepted seeds per middle space-point
    /// In dense environments many seeds may be found per middle space-point
    /// Only seeds with the highest weight will be kept if this limit is reached
    unsigned int maxSeedsPerSpM = 10;
    /// Maximum limit to number of compatible space-point used in score
    /// calculation. We increase by c1 the weight calculation for each
    /// compatible space-point until we reach compatSeedLimit
    std::size_t compatSeedLimit = 2;

    /// Increment in seed weight if the number of compatible seeds is larger
    /// than numSeedIncrement, this is used in case of high occupancy scenarios
    /// if we want to increase the weight of the seed by seedWeightIncrement
    /// when the number of compatible seeds is higher than a certain value
    float seedWeightIncrement = 0;
    float numSeedIncrement = std::numeric_limits<float>::infinity();

    // Other parameters

    /// Use deltaR between top and middle SP instead of top radius to search for
    /// compatible SPs
    bool useDeltaRorTopRadius = false;

    std::shared_ptr<IExperimentCuts2> experimentCuts;

    DerivedConfig derive() const;
  };

  struct DerivedConfig : public Config {};

  struct State {
    std::vector<std::uint32_t> topSpIndexVec;
    std::vector<float> compatibleSeedR;

    std::unordered_map<SpacePointIndex2, float> bestSeedQualityMap;
  };

  struct Options {
    bool seedConfirmation = false;
    SeedConfirmationRangeConfig seedConfRange;
    std::size_t nTopSeedConf = 0;
  };

  explicit SeedFilter2(const DerivedConfig& config,
                       std::unique_ptr<const Logger> logger = getDefaultLogger(
                           "SeedFilter2", Logging::Level::INFO));

  void filter2SpFixed(const Options& options, State& state,
                      const SpacePointContainer2& spacePoints,
                      const SpacePointColumn2<float>& rColumn,
                      SpacePointIndex2 bottomSp, SpacePointIndex2 middleSp,
                      const std::vector<SpacePointIndex2>& topSpVec,
                      const std::vector<float>& invHelixDiameterVec,
                      const std::vector<float>& impactParametersVec,
                      float zOrigin,
                      CandidatesForMiddleSp2& candidatesCollector) const;

  void filter1SpFixed(const Options& options, State& state,
                      std::vector<TripletCandidate2>& candidates,
                      std::size_t numQualitySeeds,
                      SeedContainer2& outputCollection) const;

 private:
  const Logger& logger() const { return *m_logger; }

  const DerivedConfig m_cfg;

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
