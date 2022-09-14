// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
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

template <typename external_spacepoint_t, typename platform_t = void*>
class Seedfinder {
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

 public:
  struct State {
    // bottom space point
    std::vector<InternalSpacePoint<external_spacepoint_t>*> compatBottomSP;
    std::vector<InternalSpacePoint<external_spacepoint_t>*> compatTopSP;
    // contains parameters required to calculate circle with linear equation
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;

    // create vectors here to avoid reallocation in each loop
    std::vector<InternalSpacePoint<external_spacepoint_t>*> topSpVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    std::vector<float> etaVec;
    std::vector<float> ptVec;

    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSpM;
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
  /// @param state State object that holds memory used
  /// @param outIt Output iterator for the seeds in the group
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @param rRangeSPExtent extent containing r values of all SP.
  /// @note Ranges must return pointers.
  /// @note Ranges must be separate objects for each parallel call.
  template <template <typename...> typename container_t, typename sp_range_t>
  void createSeedsForGroup(
      State& state,
      std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
      sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
      Extent rRangeSPExtent) const;

  /// @brief Compatibility method for the new-style seed finding API.
  ///
  /// This method models the old-style seeding API where we only need a
  /// container for the bottom, middle, and top space points. Also, the results
  /// are returned by value instead of inserted into an inserter.
  ///
  /// @note This method is a very simply wrapper around the more modern API.
  /// @warning The performance of the seeding code is far greater if the new
  /// API is used, and this is recommended for all new uses which do not
  /// require backwards-compatibility.
  ///
  /// @tparam sp_range_t container type for the seed point collections.
  /// @param bottomSPs group of space points to be used as innermost SP in a
  /// seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @returns a vector of seeds.
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t>> createSeedsForGroup(
      sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:
  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
};

}  // namespace Acts

#ifndef DOXYGEN
#include "Acts/Seeding/Seedfinder.ipp"
#endif
