// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointMutableData.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/IExperimentCuts.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

namespace Acts {

/// @brief Structure that holds configuration parameters for the seed filter algorithm
struct SeedFilterConfig {
  /// Allowed difference in curvature (inverted seed radii) between two
  /// compatible seeds
  float deltaInvHelixDiameter = 0.00003 * 1. / UnitConstants::mm;
  /// Minimum distance between compatible outer space-points to be considered.
  /// This is used to avoid counting space-points from the same layer
  float deltaRMin = 5. * UnitConstants::mm;
  /// Seed weight/score is increased by this value if a compatible seed has been
  /// found. This is the c1 factor in the seed score calculation (w = c1 * Nt -
  /// c2 * d0 - c3 * z0)
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
  /// calculation. We increase by c1 the weight calculation for each compatible
  /// space-point until we reach compatSeedLimit
  std::size_t compatSeedLimit = 2;

  /// Increment in seed weight if the number of compatible seeds is larger than
  /// numSeedIncrement, this is used in case of high occupancy scenarios if we
  /// want to increase the weight of the seed by seedWeightIncrement when the
  /// number of compatible seeds is higher than a certain value
  float seedWeightIncrement = 0;
  float numSeedIncrement = std::numeric_limits<float>::infinity();

  /// Seeding parameters used for quality seed confirmation

  /// Enable quality seed confirmation, this is different than default seeding
  /// confirmation because it can also be defined for different (r, z) regions
  /// of the detector (e.g. forward or central region) by SeedConfirmationRange.
  /// Seeds are classified as "high-quality" seeds and normal quality seeds.
  /// Normal quality seeds are only selected if no other "high-quality" seed
  /// has been found for that inner-middle doublet.
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
  /// Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
  /// seed confirmation. If the limit is reached we check if there is a normal
  /// quality seed to be replaced
  std::size_t maxQualitySeedsPerSpMConf =
      std::numeric_limits<std::size_t>::max();

  /// Other parameters

  /// Use deltaR between top and middle SP instead of top radius to search for
  /// compatible SPs
  bool useDeltaRorTopRadius = false;

  bool isInInternalUnits = false;
  SeedFilterConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SeedFilterConfig");
    }
    using namespace UnitLiterals;
    SeedFilterConfig config = *this;
    config.isInInternalUnits = true;
    config.deltaRMin /= 1_mm;
    config.deltaInvHelixDiameter /= 1. / 1_mm;

    return config;
  }
};

struct SeedFilterState {
  // longitudinal impact parameter as defined by bottom and middle space point
  float zOrigin = 0;
  // number of minimum top SPs in seed confirmation
  std::size_t nTopSeedConf = 0;
  // radius of bottom component of seed that is used to define the number of
  // compatible top required
  float rMaxSeedConf = std::numeric_limits<float>::max();  // UnitConstants::mm
};

/// Filter seeds at various stages with the currently
/// available information.
template <typename external_spacepoint_t>
class SeedFilter final {
 public:
  using Config = SeedFilterConfig;
  using State = SeedFilterState;

  explicit SeedFilter(
      const Config& config,
      IExperimentCuts<external_spacepoint_t>* expCuts = nullptr);
  explicit SeedFilter(
      const Config& config, std::unique_ptr<const Logger> logger,
      IExperimentCuts<external_spacepoint_t>* expCuts = nullptr);

  /// Create Seeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param bottomSP fixed bottom space point
  /// @param middleSP fixed middle space point
  /// @param topSpVec vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec vector containing 1/(2*r) values where r is the helix radius
  /// @param impactParametersVec vector containing the impact parameters
  /// @param state holds quantities used in seed filter
  /// @param candidatesCollector container for the seed candidates
  void filterSeeds_2SpFixed(
      const SpacePointMutableData& mutableData,
      const external_spacepoint_t& bottomSP,
      const external_spacepoint_t& middleSP,
      const std::vector<const external_spacepoint_t*>& topSpVec,
      const std::vector<float>& invHelixDiameterVec,
      const std::vector<float>& impactParametersVec, State& state,
      CandidatesForMiddleSp<const external_spacepoint_t>& candidatesCollector)
      const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates_collector collection of seed candidates
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename collection_t>
  void filterSeeds_1SpFixed(
      SpacePointMutableData& mutableData,
      CandidatesForMiddleSp<const external_spacepoint_t>& candidates_collector,
      collection_t& outputCollection) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates collection of seed candidates
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename collection_t>
  void filterSeeds_1SpFixed(
      SpacePointMutableData& mutableData,
      std::vector<typename CandidatesForMiddleSp<
          const external_spacepoint_t>::value_type>& candidates,
      const std::size_t numQualitySeeds, collection_t& outputCollection) const;

  const Config getConfig() const { return m_cfg; }
  const IExperimentCuts<external_spacepoint_t>* getExperimentCuts() const {
    return m_experimentCuts;
  }

 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger =
      getDefaultLogger("Filter", Logging::Level::INFO);
  const IExperimentCuts<external_spacepoint_t>* m_experimentCuts;
};

}  // namespace Acts

#include "Acts/Seeding/SeedFilter.ipp"
