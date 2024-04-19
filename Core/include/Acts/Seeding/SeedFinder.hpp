// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Neighbour.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"

#include <array>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

enum class SpacePointCandidateType : short { eBottom, eTop };

enum class DetectorMeasurementInfo : short { eDefault, eDetailed };

template <typename external_spacepoint_t, typename grid_t,
          typename platform_t = void*>
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
    std::vector<const InternalSpacePoint<external_spacepoint_t>*> topSpVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    // managing seed candidates for SpM
    CandidatesForMiddleSp<const InternalSpacePoint<external_spacepoint_t>>
        candidates_collector;

    // managing doublet candidates
    boost::container::small_vector<Acts::Neighbour<grid_t>,
                                   Acts::detail::ipow(3, grid_t::DIM)>
        bottomNeighbours;
    boost::container::small_vector<Acts::Neighbour<grid_t>,
                                   Acts::detail::ipow(3, grid_t::DIM)>
        topNeighbours;

    // Adding space point info
    Acts::SpacePointData spacePointData;
  };

  /// The only constructor. Requires a config object.
  /// @param config the configuration for the SeedFinder
  SeedFinder(const Acts::SeedFinderConfig<external_spacepoint_t>& config);
  ~SeedFinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  SeedFinder() = default;
  SeedFinder(const SeedFinder<external_spacepoint_t, grid_t, platform_t>&) =
      delete;
  SeedFinder<external_spacepoint_t, grid_t, platform_t>& operator=(
      const SeedFinder<external_spacepoint_t, grid_t, platform_t>&) = default;
  //@}

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param options frequently changing configuration (like beam position)
  /// @param state State object that holds memory used
  /// @param grid The grid with space points
  /// @param outIt Output iterator for the seeds in the group
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin.
  /// @note Ranges must return pointers.
  /// @note Ranges must be separate objects for each parallel call.
  template <typename sp_range_t>
  void createSeedsForGroup(
      const Acts::SeedFinderOptions& options, SeedingState& state,
      const grid_t& grid,
      GenericBackInserter<Seed<external_spacepoint_t>> outIt,
      const sp_range_t& bottomSPs, const std::size_t middleSPs,
      const sp_range_t& topSPs,
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
  /// @param grid The grid with space points
  /// @param bottomSPs group of space points to be used as innermost SP in a
  /// seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @returns a vector of seeds.
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t>> createSeedsForGroup(
      const Acts::SeedFinderOptions& options, const grid_t& grid,
      const sp_range_t& bottomSPs, const std::size_t middleSPs,
      const sp_range_t& topSPs) const;

 private:
  /// Iterates over dublets and tests the compatibility between them by applying
  /// a series of cuts that can be tested with only two SPs
  /// @param spacePointData object containing the spacepoint data
  /// @param options frequently changing configuration (like beam position)
  /// @param grid spacepoint grid
  /// @param otherSPsNeighbours inner or outer space points to be used in the dublet
  /// @param mediumSP space point candidate to be used as middle SP in a seed
  /// @param linCircleVec vector containing inner or outer SP parameters after reference frame transformation to the u-v space
  /// @param outVec Output object containing top or bottom SPs that are compatible with a certain middle SPs
  /// @param deltaRMinSP minimum allowed r-distance between dublet components
  /// @param deltaRMaxSP maximum allowed r-distance between dublet components
  /// @param uIP minus one over radius of middle SP
  /// @param uIP2 square of uIP
  /// @param cosPhiM ratio between middle SP x position and radius
  /// @param sinPhiM ratio between middle SP y position and radius
  template <Acts::SpacePointCandidateType candidateType, typename out_range_t>
  void getCompatibleDoublets(
      Acts::SpacePointData& spacePointData,
      const Acts::SeedFinderOptions& options, const grid_t& grid,
      boost::container::small_vector<Acts::Neighbour<grid_t>,
                                     Acts::detail::ipow(3, grid_t::DIM)>&
          otherSPsNeighbours,
      const InternalSpacePoint<external_spacepoint_t>& mediumSP,
      std::vector<LinCircle>& linCircleVec, out_range_t& outVec,
      const float deltaRMinSP, const float deltaRMaxSP, const float uIP,
      const float uIP2, const float cosPhiM, const float sinPhiM) const;

  /// Iterates over the seed candidates tests the compatibility between three
  /// SPs and calls for the seed confirmation
  /// @param spacePointData object containing the spacepoint data
  /// @param SpM space point candidate to be used as middle SP in a seed
  /// @param options frequently changing configuration (like beam position)
  /// @param seedFilterState State object that holds memory used in SeedFilter
  /// @param state State object that holds memory used
  template <Acts::DetectorMeasurementInfo detailedMeasurement>
  void filterCandidates(Acts::SpacePointData& spacePointData,
                        const InternalSpacePoint<external_spacepoint_t>& SpM,
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
