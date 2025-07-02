// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace Acts::Experimental {

class ITripletSeedCuts;

/// @brief Triplet seed filter used in the triplet seeding algorithm
///
/// Note that this algorithm is designed and tuned for cylindrical detectors and
/// uses R-Z coordinates for the space points.
class BroadTripletSeedFilter final {
 public:
  /// @brief Structure that holds configuration parameters for the seed filter algorithm
  struct Config {
    /// Allowed difference in curvature (inverted seed radii) between two
    /// compatible seeds
    float deltaInvHelixDiameter = 0.00003 * (1 / UnitConstants::mm);
    /// Minimum distance between compatible outer space-points to be considered.
    /// This is used to avoid counting space-points from the same layer
    float deltaRMin = 5 * UnitConstants::mm;
    /// Seed weight/score is increased by this value if a compatible seed has
    /// been found. This is the c1 factor in the seed score calculation (w = c1
    /// * Nt - c2 * d0 - c3 * z0)
    float compatSeedWeight = 200;
    /// The transverse impact parameters (d0) is multiplied by this factor and
    /// subtracted from weight. This is the c2 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float impactWeightFactor = 1;
    /// The logitudinal impact parameters (z0) is multiplied by this factor and
    /// subtracted from weight. This is the c3 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float zOriginWeightFactor = 1;
    /// Maximum number (minus one) of accepted seeds per middle space-point
    /// In dense environments many seeds may be found per middle space-point
    /// Only seeds with the highest weight will be kept if this limit is reached
    unsigned int maxSeedsPerSpM = 5;
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

    /// Seeding parameters used for quality seed confirmation

    /// Enable quality seed confirmation, this is different than default seeding
    /// confirmation because it can also be defined for different (r, z) regions
    /// of the detector (e.g. forward or central region) by
    /// SeedConfirmationRange. Seeds are classified as "high-quality" seeds and
    /// normal quality seeds. Normal quality seeds are only selected if no other
    /// "high-quality" seed has been found for that inner-middle doublet.
    bool seedConfirmation = false;
    /// Contains parameters for central seed confirmation
    SeedConfirmationRangeConfig centralSeedConfirmationRange;
    /// Contains parameters for forward seed confirmation
    SeedConfirmationRangeConfig forwardSeedConfirmationRange;

    /// If seedConfirmation is true we classify seeds as "high-quality" seeds.
    /// Seeds that are not confirmed as "high-quality" are only selected if no
    /// other "high-quality" seed has been found for that inner-middle doublet
    /// Maximum number of normal seeds (not classified as "high-quality" seeds)
    /// in seed confirmation
    std::size_t maxSeedsPerSpMConf = std::numeric_limits<std::size_t>::max();
    /// Maximum number of "high-quality" seeds for each inner-middle SP-dublet
    /// in seed confirmation. If the limit is reached we check if there is a
    /// normal quality seed to be replaced
    std::size_t maxQualitySeedsPerSpMConf =
        std::numeric_limits<std::size_t>::max();

    // Other parameters

    /// Use deltaR between top and middle SP instead of top radius to search for
    /// compatible SPs
    bool useDeltaRinsteadOfTopRadius = false;

    std::shared_ptr<ITripletSeedCuts> experimentCuts;
  };

  struct State {
    std::unordered_map<SpacePointIndex2, float> bestSeedQualityMap;
  };

  struct Cache {
    std::vector<std::uint32_t> topSpIndexVec;
    std::vector<float> compatibleSeedR;
  };

  explicit BroadTripletSeedFilter(const Config& config,
                                  std::unique_ptr<const Logger> logger =
                                      getDefaultLogger("BroadTripletSeedFilter",
                                                       Logging::Level::INFO));

  const Config& config() const { return m_cfg; }

  /// Create seed candidates with fixed bottom and middle space points and
  /// all compatible top space points.
  ///
  /// @param state Mutable state that is used to store intermediate results
  /// @param cache Cache object to store intermediate results
  /// @param spacePoints Container with all space points
  /// @param bottomSp Fixed bottom space point
  /// @param middleSp Fixed middle space point
  /// @param topSpVec Vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec Vector containing 1/(2*r) values where r is the
  ///                            helix radius
  /// @param impactParametersVec Vector containing the impact parameters
  /// @param zOrigin Z origin of the detector, used for z0 calculation
  /// @param candidatesCollector Container for the seed candidates
  void filter2SpFixed(State& state, Cache& cache,
                      const SpacePointContainer2& spacePoints,
                      SpacePointIndex2 bottomSp, SpacePointIndex2 middleSp,
                      std::span<const SpacePointIndex2> topSpVec,
                      std::span<const float> invHelixDiameterVec,
                      std::span<const float> impactParametersVec, float zOrigin,
                      CandidatesForMiddleSp2& candidatesCollector) const;

  /// Create final seeds for all candidates with the same middle space point
  ///
  /// @param state Mutable state that is used to store intermediate results
  /// @param spacePoints Container with all space points
  /// @param candidates Collection of seed candidates
  /// @param numQualitySeeds Number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  void filter1SpFixed(State& state, const SpacePointContainer2& spacePoints,
                      std::span<TripletCandidate2> candidates,
                      std::size_t numQualitySeeds,
                      SeedContainer2& outputCollection) const;

 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts::Experimental
