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
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/ITripletSeedFilter.hpp"
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
class BroadTripletSeedFilter final : public ITripletSeedFilter {
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
    float rMaxSeedConf;

    /// Map to store the best seed quality for each space point
    /// This is used to avoid creating seeds with lower quality than the best
    /// seed quality already found for that space point
    /// The key is the space point index, and the value is the best seed quality
    /// found for that space point
    /// @note The index is the space point index, not the seed index.
    ///       `copyFromIndex` will be used if available.
    std::unordered_map<SpacePointIndex2, float> bestSeedQualityMap;
  };

  struct Cache {
    std::vector<std::uint32_t> topSpIndexVec;
    std::vector<float> compatibleSeedR;
  };

  /// @param config Configuration parameters for the seed filter
  /// @param state Mutable state that is used to store intermediate results
  /// @param cache Cache object to store intermediate results
  /// @param logger Logger for debugging and information messages
  /// @note objects from this class depend on a per-event state and cache
  ///       and should not be used across events.
  explicit BroadTripletSeedFilter(const Config& config, State& state,
                                  Cache& cache, const Logger& logger);

  const Config& config() const { return *m_cfg; }
  State& state() const { return *m_state; }
  Cache& cache() const { return *m_cache; }
  const Logger& logger() const { return *m_logger; }

  void initialize(CandidatesForMiddleSp2& candidatesCollector) const override;

  bool sufficientTopDoublets(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp& topDoublets) const override;

  void filterTripletTopCandidates(
      const SpacePointContainer2& spacePoints,
      const DoubletsForMiddleSp::Proxy& bottomLink,
      const ConstSpacePointProxy2& spM,
      std::span<const SpacePointIndex2> topSpVec,
      std::span<const float> invHelixDiameterVec,
      std::span<const float> impactParametersVec,
      CandidatesForMiddleSp2& candidatesCollector) const override;

  void filterTripletsMiddleFixed(
      const SpacePointContainer2& spacePoints,
      std::span<TripletCandidate2> candidates, std::size_t numQualitySeeds,
      SeedContainer2& outputCollection) const override;

 private:
  const Config* m_cfg;
  State* m_state;
  Cache* m_cache;
  const Logger* m_logger;
};

}  // namespace Acts::Experimental
