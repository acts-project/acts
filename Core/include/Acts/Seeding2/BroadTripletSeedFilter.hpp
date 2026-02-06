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

namespace Acts {

/// @c ITripletSeedCuts can be used to increase or decrease seed weights
/// based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed weight cut. As the weight is stored in
/// seeds, there are two distinct methods.
class ITripletSeedCuts {
 public:
  virtual ~ITripletSeedCuts() = default;

  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  virtual float seedWeight(const ConstSpacePointProxy2& bottom,
                           const ConstSpacePointProxy2& middle,
                           const ConstSpacePointProxy2& top) const = 0;

  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  virtual bool singleSeedCut(float weight, const ConstSpacePointProxy2& bottom,
                             const ConstSpacePointProxy2& middle,
                             const ConstSpacePointProxy2& top) const = 0;

  /// @param seedCandidates contains collection of seed candidates created for one middle
  /// space point in a std::tuple format
  virtual void cutPerMiddleSp(
      std::span<TripletCandidate2>& seedCandidates) const = 0;
};

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
    /// Number of seeds required before `seedWeightIncrement` is applied
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
    std::uint32_t maxSeedsPerSpMConf = 5;
    /// Maximum number of "high-quality" seeds for each inner-middle SP-dublet
    /// in seed confirmation. If the limit is reached we check if there is a
    /// normal quality seed to be replaced
    std::uint32_t maxQualitySeedsPerSpMConf = 5;

    // Other parameters

    /// Use deltaR between top and middle SP instead of top radius to search for
    /// compatible SPs
    bool useDeltaRinsteadOfTopRadius = false;

    /// Custom cuts interface for experiments
    std::shared_ptr<ITripletSeedCuts> experimentCuts;
  };

  /// State of the filter that is communicated between different stages
  struct State {
    /// Collector for triplet candidates associated with middle space points
    CandidatesForMiddleSp2 candidatesCollector;

    /// Maximum radius for seed confirmation
    float rMaxSeedConf{};

    /// Map to store the best seed quality for each space point
    /// This is used to avoid creating seeds with lower quality than the best
    /// seed quality already found for that space point
    /// The key is the space point index, and the value is the best seed quality
    /// found for that space point
    /// @note The index is the space point index, not the seed index.
    ///       `copyFromIndex` will be used if available.
    std::unordered_map<SpacePointIndex2, float> bestSeedQualityMap;
  };

  /// Cache for intermediate results to avoid reallocations. No information is
  /// carried over between different stages.
  struct Cache {
    /// Cache for top space point indices during compatibility search
    std::vector<std::uint32_t> topSpIndexVec;
    /// Cache for compatible seed radii during score calculation
    std::vector<float> compatibleSeedR;
    /// Cache for sorted triplet candidates during selection
    std::vector<TripletCandidate2> sortedCandidates;
  };

  /// @param config Configuration parameters for the seed filter
  /// @param state Mutable state that is used to store intermediate results
  /// @param cache Cache object to store intermediate results
  /// @param logger Logger for debugging and information messages
  /// @note objects from this class depend on a per-event state and cache
  ///       and should not be used across events.
  explicit BroadTripletSeedFilter(const Config& config, State& state,
                                  Cache& cache, const Logger& logger);

  /// @param spacePoints Container of space points
  /// @param spM Middle space point proxy
  /// @param topDoublets Collection of top doublets for the middle space point
  /// @return true if sufficient top doublets are found
  bool sufficientTopDoublets(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp& topDoublets) const override;

  void filterTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomLink,
      const TripletTopCandidates& tripletTopCandidates) const override;

  void filterTripletsMiddleFixed(
      const SpacePointContainer2& spacePoints,
      SeedContainer2& outputCollection) const override;

 private:
  /// Configuration parameters for the seed filter algorithm
  const Config* m_cfg{};
  /// Mutable state for intermediate results between filter stages
  State* m_state{};
  /// Cache for temporary data to avoid reallocations
  Cache* m_cache{};
  /// Logger for debugging and information messages
  const Logger* m_logger{};

  const Config& config() const { return *m_cfg; }
  State& state() const { return *m_state; }
  Cache& cache() const { return *m_cache; }
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
