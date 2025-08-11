// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts::Experimental {

/// @brief Triplet seeding algorithm front-end
///
/// This class implements the triplet seeding algorithm, which is typical
/// procedure to find track seeds using space points in a cylindrical detector.
/// It is designed to be fast, flexible, and configurable.
///
/// The algorithm works by first finding compatible doublets of space points,
/// two space points that can be connected by a track coming from the
/// interaction region, and then forming triplets by combinding these
/// doublets at a common middle space point. The triplets are then filtered
/// using a seed filter to produce a set of track seeds.
///
/// Note that this algorithm is designed and tuned for cylindrical detectors and
/// uses R-Z coordinates for the space points.
class BroadTripletSeedFinder {
 public:
  struct Cache {
    DoubletsForMiddleSp bottomDoublets;
    DoubletsForMiddleSp topDoublets;

    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> sortedBottoms;
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> sortedTops;

    TripletTopCandidates tripletTopCandidates;

    BroadTripletSeedFilter::Cache filter;

    CandidatesForMiddleSp2 candidatesCollector;
    std::vector<TripletCandidate2> sortedCandidates;
  };

  struct State {
    BroadTripletSeedFilter::State filter;
  };

  explicit BroadTripletSeedFinder(std::unique_ptr<const Logger> logger =
                                      getDefaultLogger("BroadTripletSeedFinder",
                                                       Logging::Level::INFO));

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param bottomFinder Finder for bottom doublets
  /// @param topFinder Finder for top doublets
  /// @param tripletFinder Finder for triplet space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param spacePoints Space point container
  /// @param bottomSps Subset of space points to be used as innermost SP in a seed
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param topSps Subset of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroup(State& state, Cache& cache,
                            const DoubletSeedFinder& bottomFinder,
                            const DoubletSeedFinder& topFinder,
                            const TripletSeedFinder& tripletFinder,
                            const BroadTripletSeedFilter& filter,
                            const SpacePointContainer2& spacePoints,
                            SpacePointContainer2::ConstSubset& bottomSps,
                            const ConstSpacePointProxy2& middleSp,
                            SpacePointContainer2::ConstSubset& topSps,
                            SeedContainer2& outputSeeds) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param bottomFinder Finder for bottom doublets
  /// @param topFinder Finder for top doublets
  /// @param tripletFinder Finder for triplet space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param spacePoints Space point container
  /// @param bottomSpGroups Groups of space points to be used as innermost SP in a seed
  /// @param middleSpGroup Group of space points to be used as middle SP in a seed
  /// @param topSpGroups Groups of space points to be used as outermost SP in a seed
  /// @param radiusRangeForMiddle Range of radii for the middle space points
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroups(
      State& state, Cache& cache, const DoubletSeedFinder& bottomFinder,
      const DoubletSeedFinder& topFinder,
      const TripletSeedFinder& tripletFinder,
      const BroadTripletSeedFilter& filter,
      const SpacePointContainer2& spacePoints,
      const std::span<SpacePointContainer2::ConstRange>& bottomSpGroups,
      const SpacePointContainer2::ConstRange& middleSpGroup,
      const std::span<SpacePointContainer2::ConstRange>& topSpGroups,
      const std::pair<float, float>& radiusRangeForMiddle,
      SeedContainer2& outputSeeds) const;

 private:
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
