// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <mutex>
#include <queue>
#include <vector>

#include "Acts/Seeding/IExperimentCuts.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"

namespace Acts {
struct SeedFilterConfig {
  // the allowed delta between two inverted seed radii for them to be considered
  // compatible.
  float deltaInvHelixDiameter = 0.00003;
  // the impact parameters (d0) is multiplied by this factor and subtracted from
  // weight
  float impactWeightFactor = 1.;
  // seed weight increased by this value if a compatible seed has been found.
  float compatSeedWeight = 200.;
  // minimum distance between compatible seeds to be considered for weight boost
  float deltaRMin = 5.;
  // in dense environments many seeds may be found per middle space point.
  // only seeds with the highest weight will be kept if this limit is reached.
  unsigned int maxSeedsPerSpM = 10;
  // how often do you want to increase the weight of a seed for finding a
  // compatible seed?
  size_t compatSeedLimit = 2;
  // Tool to apply experiment specific cuts on collected middle space points
};

/// @class Filter seeds at various stages with the currently
/// available information.
template <typename external_spacepoint_t>
class SeedFilter {
 public:
  SeedFilter(SeedFilterConfig config,
             IExperimentCuts<external_spacepoint_t>* expCuts = 0);

  SeedFilter() = delete;
  virtual ~SeedFilter() = default;

  /// Create InternalSeeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param bottomSP fixed bottom space point
  /// @param middleSP fixed middle space point
  /// @param topSpVec vector containing all space points that may be compatible
  /// with both bottom and middle space point
  /// @param origin on the z axis as defined by bottom and middle space point
  /// @return vector of pairs containing seed weight and seed for all valid
  /// created seeds
  virtual std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
  filterSeeds_2SpFixed(
      const InternalSpacePoint<external_spacepoint_t>& bottomSP,
      const InternalSpacePoint<external_spacepoint_t>& middleSP,
      std::vector<const InternalSpacePoint<external_spacepoint_t>*>& topSpVec,
      std::vector<float>& invHelixDiameterVec,
      std::vector<float>& impactParametersVec, float zOrigin) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param seedsPerSpM vector of pairs containing weight and seed for all
  /// for all seeds with the same middle space point
  /// @return vector of all InternalSeeds that not filtered out
  virtual void filterSeeds_1SpFixed(
      std::vector<std::pair<
          float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>&
          seedsPerSpM,
      std::vector<Seed<external_spacepoint_t>>& outVec) const;
  const SeedFilterConfig getSeedFilterConfig() const { return m_cfg; }
  const IExperimentCuts<external_spacepoint_t>* getExperimentCuts() const {
    return m_experimentCuts;
  }

 private:
  const SeedFilterConfig m_cfg;
  const IExperimentCuts<external_spacepoint_t>* m_experimentCuts;
};
}  // namespace Acts
#include "Acts/Seeding/SeedFilter.ipp"
