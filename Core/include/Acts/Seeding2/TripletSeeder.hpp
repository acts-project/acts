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
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/ITripletSeedFilter.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts {

/// Full triplet seeder which depends on a doublet and triplet seed finder, and
/// a triplet seed filter.
class TripletSeeder {
 public:
  /// Cache for storing intermediate results during triplet seeding to avoid
  /// reallocation.
  struct Cache {
    /// Cache for bottom doublets associated with middle space points
    DoubletsForMiddleSp bottomDoublets;
    /// Cache for top doublets associated with middle space points
    DoubletsForMiddleSp topDoublets;

    /// Sorted container of bottom doublet indices with cotangent theta values
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> sortedBottoms;
    /// Sorted container of top doublet indices with cotangent theta values
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> sortedTops;

    /// Cache for triplet top candidates during seed formation
    TripletTopCandidates tripletTopCandidates;
  };

  /// Construct a TripletSeeder with optional logger.
  /// @param logger Logger instance for debug output (defaults to INFO level)
  explicit TripletSeeder(std::unique_ptr<const Logger> logger =
                             getDefaultLogger("TripletSeeder",
                                              Logging::Level::INFO));

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
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
  void createSeedsFromGroup(Cache& cache, const DoubletSeedFinder& bottomFinder,
                            const DoubletSeedFinder& topFinder,
                            const TripletSeedFinder& tripletFinder,
                            const ITripletSeedFilter& filter,
                            const SpacePointContainer2& spacePoints,
                            SpacePointContainer2::ConstSubset& bottomSps,
                            const ConstSpacePointProxy2& middleSp,
                            SpacePointContainer2::ConstSubset& topSps,
                            SeedContainer2& outputSeeds) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
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
      Cache& cache, const DoubletSeedFinder& bottomFinder,
      const DoubletSeedFinder& topFinder,
      const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
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

}  // namespace Acts
