// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/TripletSeeder.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <span>

#include <Eigen/Dense>

namespace Acts {

namespace {

/// Advance a space point group past the candidates that can no longer pair
/// with the middle space point, or with any of the ones after it. Groups are
/// sorted by radius and the middle space points are swept in ascending radius,
/// so the cut point only ever moves forward and the group is never rewound.
///
/// This is the incremental form of the binary search that primes the groups in
/// `createSeedsFromGroups`; `isInRange` has to spell the radius cut exactly the
/// way `DoubletSeedFinder` does, so that the two agree bit for bit.
template <typename SpacePointGroup, typename Predicate>
std::uint32_t advanceGroup(SpacePointGroup& group, Predicate isInRange) {
  std::uint32_t offset = 0;
  for (ConstSpacePointProxy candidateSp : group) {
    if (isInRange(candidateSp)) {
      break;
    }
    ++offset;
  }
  group = group.subrange(offset);
  return offset;
}

/// Move `end` forward past the candidates in the middle space point's radial
/// range; `isInRange` must spell the radius cut as `DoubletSeedFinder` does.
template <typename SpacePointGroup, typename Predicate>
std::uint32_t advanceGroupEnd(const SpacePointGroup& group, std::uint32_t end,
                              Predicate isInRange) {
  const auto size = static_cast<std::uint32_t>(group.size());
  while (end < size && isInRange(group[end])) {
    ++end;
  }
  return end;
}

/// The columns kept longer than the candidates by this many values, which
/// the doublet finders read to make their flags for whole vectors: zeros, or
/// candidates of an earlier batch
constexpr std::size_t kColumnPadding = 16;

/// Make room in the columns for a number of candidates
void reserveColumns(TripletSeeder::Cache::CandidateColumns& columns,
                    std::size_t size) {
  if (columns.z.size() < size + kColumnPadding) {
    columns.z.resize(size + kColumnPadding);
    columns.r.resize(size + kColumnPadding);
    columns.indices.resize(size + kColumnPadding);
  }
}

/// The seeds of one middle space point. The candidates are what the doublet
/// finders take: they make the radial cut, so the candidates may be more than
/// the middle space point's own radial range holds.
template <typename CandidateSps>
void createSeedsFromGroupsImpl(
    const Logger& logger, TripletSeeder::Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
    const SpacePointContainer& spacePoints, const CandidateSps& bottomSps,
    const ConstSpacePointProxy& middleSp, const CandidateSps& topSps,
    SeedContainer& outputSeeds) {
  MiddleSpInfo middleSpInfo = DoubletSeedFinder::computeMiddleSpInfo(middleSp);

  // create middle-top doublets
  cache.topDoublets.clear();
  topFinder.createDoublets(middleSp, middleSpInfo, topSps, cache.topDoublets);

  // no top SP found -> cannot form any triplet
  if (cache.topDoublets.empty()) {
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  if (!filter.sufficientTopDoublets(spacePoints, middleSp, cache.topDoublets)) {
    return;
  }

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  bottomFinder.createDoublets(middleSp, middleSpInfo, bottomSps,
                              cache.bottomDoublets);

  // no bottom SP found -> cannot form any triplet
  if (cache.bottomDoublets.empty()) {
    ACTS_VERBOSE("No compatible Bottoms, returning");
    return;
  }

  ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size() << " bottoms and "
                              << cache.topDoublets.size()
                              << " tops for middle candidate indexed "
                              << middleSp.index());

  // the triplet top candidates of each bottom doublet go to the filter
  const auto filterCandidates =
      [&](const DoubletsForMiddleSp::Proxy& bottomDoublet,
          const TripletTopCandidates& candidates) {
        filter.filterTripletTopCandidates(spacePoints, middleSp, bottomDoublet,
                                          candidates);
      };
  TripletSeedFinder::CandidatesCallback filterCallback;
  filterCallback.connect(filterCandidates);

  // combine doublets to triplets
  if (tripletFinder.config().sortedByCotTheta) {
    cache.bottomDoublets.sortByCotTheta({0, cache.bottomDoublets.size()},
                                        cache.sortedBottoms, cache.sortStorage);
    cache.topDoublets.sortByCotTheta({0, cache.topDoublets.size()},
                                     cache.sortedTops, cache.sortStorage);

    tripletFinder.createTripletTopCandidates(
        spacePoints, middleSp, cache.bottomDoublets.subset(cache.sortedBottoms),
        cache.topDoublets.subset(cache.sortedTops), cache.tripletTopCandidates,
        filterCallback);
  } else {
    tripletFinder.createTripletTopCandidates(
        spacePoints, middleSp, cache.bottomDoublets.range(),
        cache.topDoublets.range(), cache.tripletTopCandidates, filterCallback);
  }

  filter.filterTripletsMiddleFixed(spacePoints, outputSeeds);
}

}  // namespace

TripletSeeder::TripletSeeder(std::unique_ptr<const Logger> logger_)
    : m_logger(std::move(logger_)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument("TripletSeeder: logger cannot be null");
  }
}

void TripletSeeder::createSeedsFromGroup(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer& spacePoints,
    SpacePointContainer::ConstSubset bottomSps,
    const ConstSpacePointProxy& middleSp,
    SpacePointContainer::ConstSubset topSps, SeedContainer& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");

  createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                            tripletFinder, filter, spacePoints, bottomSps,
                            middleSp, topSps, outputSeeds);
}

void TripletSeeder::createSeedsFromGroups(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer& spacePoints,
    std::span<const SpacePointContainer::ConstRange> bottomSpGroups,
    const SpacePointContainer::ConstRange& middleSpGroup,
    std::span<const SpacePointContainer::ConstRange> topSpGroups,
    const std::pair<float, float>& radiusRangeForMiddle,
    SeedContainer& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");
  const bool spacePointsSortedByRadius =
      bottomFinder.config().spacePointsSortedByRadius;
  const float bottomDeltaRMin = bottomFinder.config().deltaRMin;
  const float bottomDeltaRMax = bottomFinder.config().deltaRMax;
  const float topDeltaRMin = topFinder.config().deltaRMin;
  const float topDeltaRMax = topFinder.config().deltaRMax;

  if (middleSpGroup.empty()) {
    return;
  }

  // the cursors are ours, so the caller's groups stay untouched and can be
  // handed to another call unchanged
  cache.bottomSpGroups.assign(bottomSpGroups.begin(), bottomSpGroups.end());
  cache.topSpGroups.assign(topSpGroups.begin(), topSpGroups.end());
  cache.bottomSpGroupEnds.assign(bottomSpGroups.size(), 0);
  cache.topSpGroupEnds.assign(topSpGroups.size(), 0);
  // what the doublet finders get: the whole groups, or, sorted by radius,
  // their part in the radial range of a few middle space points
  cache.bottomCandidateSpGroups = cache.bottomSpGroups;
  cache.topCandidateSpGroups = cache.topSpGroups;

  // each group advanced to the first space point in the radial range of a
  // middle space point and cut at the end of the range of a later one
  const auto cutGroups = [](auto& spGroups, auto& ends, auto& candidateSpGroups,
                            auto reachesWindow, auto withinWindow) {
    for (std::size_t g = 0; g < spGroups.size(); ++g) {
      auto& spGroup = spGroups[g];
      const std::uint32_t offset = advanceGroup(spGroup, reachesWindow);
      // the start may have passed the old end of an empty window
      ends[g] = advanceGroupEnd(spGroup, ends[g] - std::min(ends[g], offset),
                                withinWindow);
      candidateSpGroups[g] = spGroup.subrange(0, ends[g]);
    }
  };
  // the candidates of all groups copied into columns, one group after the
  // other: the doublet finders then go through them in one piece, and read
  // them from a few hundred bytes for each middle space point of the batch
  const auto zrColumn = spacePoints.zrColumn().data();
  const auto copyCandidates =
      [&](const std::vector<SpacePointContainer::ConstRange>& candidateSpGroups,
          Cache::CandidateColumns& columns) {
        std::size_t next = 0;
        for (const SpacePointContainer::ConstRange& spGroup :
             candidateSpGroups) {
          const SpacePointIndex first = spGroup.range().first;
          const std::size_t count = spGroup.size();
          reserveColumns(columns, next + count);
          float* const z = columns.z.data() + next;
          float* const r = columns.r.data() + next;
          SpacePointIndex* const indices = columns.indices.data() + next;
          for (std::size_t k = 0; k < count; ++k) {
            z[k] = zrColumn[first + k][0];
            r[k] = zrColumn[first + k][1];
            indices[k] = static_cast<SpacePointIndex>(first + k);
          }
          next += count;
        }
        columns.count = next;
      };

  // the middle space points that share one cut of the groups: their radial
  // ranges nearly coincide, and the doublet finders cut each one's exactly
  constexpr std::uint32_t middleSpsPerCut = 8;
  std::uint32_t middleSpsUntilCut = 0;

  if (spacePointsSortedByRadius) {
    // Initialize initial offsets for bottom and top space points with binary
    // search. This requires at least one middle space point to be present which
    // is already checked above.
    const ConstSpacePointProxy firstMiddleSp = middleSpGroup.front();
    const float firstMiddleSpR = firstMiddleSp.zr()[1];

    for (auto& bottomSpGroup : cache.bottomSpGroups) {
      // Find the first bottom space point that is within the deltaRMax of the
      // first middle space point.
      const auto low = std::ranges::lower_bound(
          bottomSpGroup, firstMiddleSpR - bottomDeltaRMax, {},
          [&](const ConstSpacePointProxy& sp) { return sp.zr()[1]; });
      bottomSpGroup = bottomSpGroup.subrange(low - bottomSpGroup.begin());
    }

    for (auto& topSpGroup : cache.topSpGroups) {
      // Find the first top space point that is within the deltaRMin of the
      // first middle space point.
      const auto low = std::ranges::lower_bound(
          topSpGroup, firstMiddleSpR + topDeltaRMin, {},
          [&](const ConstSpacePointProxy& sp) { return sp.zr()[1]; });
      topSpGroup = topSpGroup.subrange(low - topSpGroup.begin());
    }
  }

  for (std::uint32_t m = 0; m < middleSpGroup.size(); ++m) {
    const ConstSpacePointProxy spM = middleSpGroup[m];
    const float rM = spM.zr()[1];

    if (spacePointsSortedByRadius) {
      // check if spM is outside our radial region of interest
      if (rM < radiusRangeForMiddle.first) {
        continue;
      }
      if (rM > radiusRangeForMiddle.second) {
        // break because SPs are sorted in r
        break;
      }

      if (middleSpsUntilCut == 0) {
        middleSpsUntilCut = middleSpsPerCut;
        const std::uint32_t last =
            std::min<std::uint32_t>(m + middleSpsPerCut, middleSpGroup.size()) -
            1;
        const float rLast = middleSpGroup[last].zr()[1];
        cutGroups(
            cache.topSpGroups, cache.topSpGroupEnds, cache.topCandidateSpGroups,
            [&](const ConstSpacePointProxy& sp) {
              return sp.zr()[1] - rM >= topDeltaRMin;
            },
            [&](const ConstSpacePointProxy& sp) {
              return sp.zr()[1] - rLast <= topDeltaRMax;
            });
        cutGroups(
            cache.bottomSpGroups, cache.bottomSpGroupEnds,
            cache.bottomCandidateSpGroups,
            [&](const ConstSpacePointProxy& sp) {
              return rM - sp.zr()[1] <= bottomDeltaRMax;
            },
            [&](const ConstSpacePointProxy& sp) {
              return rLast - sp.zr()[1] >= bottomDeltaRMin;
            });
        copyCandidates(cache.topCandidateSpGroups, cache.topCandidates);
        copyCandidates(cache.bottomCandidateSpGroups, cache.bottomCandidates);
      }
      --middleSpsUntilCut;

      createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                                tripletFinder, filter, spacePoints,
                                cache.bottomCandidates, spM,
                                cache.topCandidates, outputSeeds);
      continue;
    }

    createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                              tripletFinder, filter, spacePoints,
                              std::span<const SpacePointContainer::ConstRange>(
                                  cache.bottomCandidateSpGroups),
                              spM,
                              std::span<const SpacePointContainer::ConstRange>(
                                  cache.topCandidateSpGroups),
                              outputSeeds);
  }
}

}  // namespace Acts
