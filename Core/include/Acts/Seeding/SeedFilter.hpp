// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/IExperimentCuts.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"

#include <memory>
#include <mutex>
#include <queue>
#include <vector>

namespace Acts {
/// Filter seeds at various stages with the currently
/// available information.

class SeedFilter {
 public:
  SeedFilter(SeedFilterConfig config, IExperimentCuts* expCuts = 0);

  SeedFilter() = delete;
  virtual ~SeedFilter() = default;

  /// Create InternalSeeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param bottomSP fixed bottom space point
  /// @param middleSP fixed middle space point
  /// @param topSpVec vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec vector containing 1/(2*r) values where r is the helix radius
  /// @param impactParametersVec vector containing the impact parameters
  /// @param zOrigin on the z axis as defined by bottom and middle space point
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param numSeeds number of seeds that did not pass the quality confirmation but were still accepted, if quality confirmation is not used this is the total number of seeds
  /// @param outCont Output container for the seeds
  virtual void filterSeeds_2SpFixed(
      Acts::SpacePoint& bottomSP, Acts::SpacePoint& middleSP,
      std::vector<Acts::SpacePoint*>& topSpVec,
      std::vector<float>& invHelixDiameterVec,
      std::vector<float>& impactParametersVec, float zOrigin,
      int& numQualitySeeds, int& numSeeds,
      std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>&
          outCont) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param seedsPerSpM vector of pairs containing weight and seed for all
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param outIt Output iterator for the seeds
  /// for all seeds with the same middle space point
  virtual void filterSeeds_1SpFixed(
      std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>&
          seedsPerSpM,
      int& numQualitySeeds,
      std::back_insert_iterator<std::vector<InternalSeed>> outIt) const;

  /// Check if there is a lower quality seed that can be replaced
  /// @param bottomSP fixed bottom space point
  /// @param middleSP fixed middle space point
  /// @param topSp fixed top space point
  /// @param zOrigin on the z axis as defined by bottom and middle space point
  /// @param isQualitySeed information whether the seed is quality confirmed or not
  /// @param weight weight of the seed
  /// @param outCont container for the seeds
  virtual void checkReplaceSeeds(
      Acts::SpacePoint& bottomSP, Acts::SpacePoint& middleSP,
      Acts::SpacePoint& topSp, float zOrigin, bool isQualitySeed, float weight,
      std::vector<std::pair<float, std::unique_ptr<const InternalSeed>>>&
          outCont) const;

  const SeedFilterConfig getSeedFilterConfig() const { return m_cfg; }
  const IExperimentCuts* getExperimentCuts() const { return m_experimentCuts; }

 private:
  const SeedFilterConfig m_cfg;
  const IExperimentCuts* m_experimentCuts;
};
}  // namespace Acts
