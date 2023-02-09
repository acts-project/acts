// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"

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
class SeedFinder {
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

 public:
  struct SeedingState {
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

    // managing seed candidates for SpM
    CandidatesForMiddleSp<InternalSpacePoint<external_spacepoint_t>>
        candidates_collector;
  };

  /// The only constructor. Requires a config object.
  /// @param config the configuration for the SeedFinder
  SeedFinder(const Acts::SeedFinderConfig<external_spacepoint_t>& config);
  ~SeedFinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  SeedFinder() = default;
  SeedFinder(const SeedFinder<external_spacepoint_t, platform_t>&) = delete;
  SeedFinder<external_spacepoint_t, platform_t>& operator=(
      const SeedFinder<external_spacepoint_t, platform_t>&) = default;
  //@}

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param options frequently changing configuration (like beam position)
  /// @param state State object that holds memory used
  /// @param outIt Output iterator for the seeds in the group
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin.
  /// @note Ranges must return pointers.
  /// @note Ranges must be separate objects for each parallel call.
  template <template <typename...> typename container_t, typename sp_range_t>
  void createSeedsForGroup(
      const Acts::SeedFinderOptions& options, SeedingState& state,
      std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
      sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
      const Acts::Range1D<float>& rMiddleSPRange) const;

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
  /// @param options frequently changing configuration (like beam position)
  /// @param bottomSPs group of space points to be used as innermost SP in a
  /// seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @returns a vector of seeds.
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t>> createSeedsForGroup(
      const Acts::SeedFinderOptions& options, sp_range_t bottomSPs,
      sp_range_t middleSPs, sp_range_t topSPs) const;

 private:
  template <typename sp_range_t, typename out_range_t>
  void getCompatibleDoublets(
      const Acts::SeedFinderOptions& options, sp_range_t& otherSPs,
      const InternalSpacePoint<external_spacepoint_t>& mediumSP,
      out_range_t& outVec, std::vector<LinCircle>& linCircleVec,
      const float& deltaRMinSP, const float& deltaRMaxSP, bool isBottom) const;

  void filterCandidates(InternalSpacePoint<external_spacepoint_t>& SpM,
                        const Acts::SeedFinderOptions& options,
                        SeedFilterState& seedFilterState,
                        SeedingState& state) const;

 private:
  Acts::SeedFinderConfig<external_spacepoint_t> m_config;
};

}  // namespace Acts

#ifndef DOXYGEN
#include "Acts/Seeding/SeedFinder.ipp"
#endif
