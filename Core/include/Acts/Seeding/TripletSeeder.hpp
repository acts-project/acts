// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/ITripletSeedFilter.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
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
    /// Storage for the sorts of the bottom and the top doublets
    DoubletsForMiddleSp::SortStorage sortStorage;

    /// Cache for triplet top candidates during seed formation
    TripletTopCandidates tripletTopCandidates;

    /// Cursors into the bottom space point groups, advanced as the middle
    /// space points are swept in ascending radius
    std::vector<SpacePointContainer::ConstRange> bottomSpGroups;
    /// Cursors into the top space point groups, advanced the same way
    std::vector<SpacePointContainer::ConstRange> topSpGroups;
    /// End of each bottom group's radial range, as a count from its cursor
    std::vector<std::uint32_t> bottomSpGroupEnds;
    /// End of each top group's radial range, as a count from its cursor
    std::vector<std::uint32_t> topSpGroupEnds;
    /// Each bottom group as passed to the doublet finder, cut to the radial
    /// range of the current batch of middle space points
    std::vector<SpacePointContainer::ConstRange> bottomCandidateSpGroups;
    /// Each top group as passed to the doublet finder, cut to the radial
    /// range of the current batch of middle space points
    std::vector<SpacePointContainer::ConstRange> topCandidateSpGroups;

    /// Columns of candidates copied out of the space point container
    using CandidateColumns = DoubletSeedFinder::CandidateColumns;
    /// The bottom candidates of the current batch of middle space points, of
    /// all groups one after the other
    CandidateColumns bottomCandidates;
    /// The top candidates of the current batch of middle space points
    CandidateColumns topCandidates;

    /// Number of parts the middle space points of a batch are split in by z
    static constexpr std::size_t kSubBatches = 4;
    /// The bottom candidates that the middle space points of each part of
    /// the batch can reach in z
    std::array<CandidateColumns, kSubBatches> bottomSubCandidates;
    /// The top candidates that the middle space points of each part of the
    /// batch can reach in z
    std::array<CandidateColumns, kSubBatches> topSubCandidates;
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
  /// @param spacePoints space point container
  /// @param bottomSps Subset of space points to be used as innermost SP in a seed
  /// @param middleSp space point candidate to be used as middle SP in a seed
  /// @param topSps Subset of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroup(Cache& cache, const DoubletSeedFinder& bottomFinder,
                            const DoubletSeedFinder& topFinder,
                            const TripletSeedFinder& tripletFinder,
                            const ITripletSeedFilter& filter,
                            const SpacePointContainer& spacePoints,
                            SpacePointContainer::ConstSubset bottomSps,
                            const ConstSpacePointProxy& middleSp,
                            SpacePointContainer::ConstSubset topSps,
                            SeedContainer& outputSeeds) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param cache Cache object to store intermediate results
  /// @param bottomFinder Finder for bottom doublets
  /// @param topFinder Finder for top doublets
  /// @param tripletFinder Finder for triplet space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param spacePoints space point container
  /// @param bottomSpGroups Groups of space points to be used as innermost SP in a seed
  /// @param middleSpGroup Group of space points to be used as middle SP in a seed
  /// @param topSpGroups Groups of space points to be used as outermost SP in a seed
  /// @param radiusRangeForMiddle Range of radii for the middle space points
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroups(
      Cache& cache, const DoubletSeedFinder& bottomFinder,
      const DoubletSeedFinder& topFinder,
      const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
      const SpacePointContainer& spacePoints,
      std::span<const SpacePointContainer::ConstRange> bottomSpGroups,
      const SpacePointContainer::ConstRange& middleSpGroup,
      std::span<const SpacePointContainer::ConstRange> topSpGroups,
      const std::pair<float, float>& radiusRangeForMiddle,
      SeedContainer& outputSeeds) const;

 private:
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
