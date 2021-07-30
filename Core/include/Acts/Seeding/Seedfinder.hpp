// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"

#include <array>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
struct LinCircle {
  float Zo;
  float cotTheta;
  float iDeltaR;
  float Er;
  float U;
  float V;
};
template <typename external_spacepoint_t, typename platform_t = void*>
class Seedfinder {
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

 public:
  struct State {
    // bottom space point
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>
        compatBottomSP;
    std::vector<const InternalSpacePoint<external_spacepoint_t>*> compatTopSP;
    // contains parameters required to calculate circle with linear equation
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;

    // create vectors here to avoid reallocation in each loop
    std::vector<const InternalSpacePoint<external_spacepoint_t>*> topSpVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSpM;

    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        sameTrackSeeds;
  };

  /// The only constructor. Requires a config object.
  /// @param config the configuration for the Seedfinder
  Seedfinder(Acts::SeedfinderConfig<external_spacepoint_t> config);
  ~Seedfinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  Seedfinder() = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t, platform_t>&) = delete;
  Seedfinder<external_spacepoint_t, platform_t>& operator=(
      const Seedfinder<external_spacepoint_t, platform_t>&) = delete;
  //@}

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @note Ranges must return pointers.
  /// @note Ranges must be separate objects for each parallel call.
  /// @return vector in which all found seeds for this group are stored.
  template <typename sp_range_t>
  void createSeedsForGroup(
      State& state,
      std::back_insert_iterator<std::vector<Seed<external_spacepoint_t>>> outIt,
      sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:
  void transformCoordinates(
      std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
      const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
      std::vector<LinCircle>& linCircleVec) const;

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
};

}  // namespace Acts

#include "Acts/Seeding/Seedfinder.ipp"
