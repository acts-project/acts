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
void advanceGroup(SpacePointGroup& group, Predicate isInRange) {
  std::uint32_t offset = 0;
  for (ConstSpacePointProxy candidateSp : group) {
    if (isInRange(candidateSp)) {
      break;
    }
    ++offset;
  }
  group = group.subrange(offset);
}

template <typename DoubletCollections>
void createAndFilterTriplets(TripletSeeder::Cache& cache,
                             const TripletSeedFinder& tripletFinder,
                             const ITripletSeedFilter& filter,
                             const SpacePointContainer& spacePoints,
                             DoubletCollections bottomDoublets,
                             const ConstSpacePointProxy& spM,
                             DoubletCollections topDoublets) {
  for (auto bottomDoublet : bottomDoublets) {
    if (topDoublets.empty()) {
      break;
    }

    cache.tripletTopCandidates.clear();
    topDoublets = tripletFinder.createTripletTopCandidates(
        spacePoints, spM, bottomDoublet, topDoublets,
        cache.tripletTopCandidates);

    filter.filterTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      cache.tripletTopCandidates);
  }
}

template <typename SpacePointCollections>
void createSeedsFromGroupsImpl(
    const Logger& logger, TripletSeeder::Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
    const SpacePointContainer& spacePoints,
    SpacePointCollections& bottomSpGroups, const ConstSpacePointProxy& middleSp,
    SpacePointCollections& topSpGroups, bool spacePointsSortedByRadius,
    float bottomDeltaRMax, float topDeltaRMin, SeedContainer& outputSeeds) {
  MiddleSpInfo middleSpInfo = DoubletSeedFinder::computeMiddleSpInfo(middleSp);
  const float rM = middleSp.zr()[1];

  // create middle-top doublets
  cache.topDoublets.clear();
  for (auto& topSpGroup : topSpGroups) {
    if (spacePointsSortedByRadius) {
      // if r-distance is too small, try next SP in bin
      advanceGroup(topSpGroup, [&](const ConstSpacePointProxy& candidateSp) {
        return candidateSp.zr()[1] - rM >= topDeltaRMin;
      });
    }
    topFinder.createDoublets(middleSp, middleSpInfo, topSpGroup,
                             cache.topDoublets);
  }

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
  for (auto& bottomSpGroup : bottomSpGroups) {
    if (spacePointsSortedByRadius) {
      // if r-distance is too big, try next SP in bin
      advanceGroup(bottomSpGroup, [&](const ConstSpacePointProxy& candidateSp) {
        return rM - candidateSp.zr()[1] <= bottomDeltaRMax;
      });
    }
    bottomFinder.createDoublets(middleSp, middleSpInfo, bottomSpGroup,
                                cache.bottomDoublets);
  }

  // no bottom SP found -> cannot form any triplet
  if (cache.bottomDoublets.empty()) {
    ACTS_VERBOSE("No compatible Bottoms, returning");
    return;
  }

  ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size() << " bottoms and "
                              << cache.topDoublets.size()
                              << " tops for middle candidate indexed "
                              << middleSp.index());

  // combine doublets to triplets
  if (tripletFinder.config().sortedByCotTheta) {
    cache.bottomDoublets.sortByCotTheta({0, cache.bottomDoublets.size()},
                                        cache.sortedBottoms);
    cache.topDoublets.sortByCotTheta({0, cache.topDoublets.size()},
                                     cache.sortedTops);

    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.subset(cache.sortedBottoms),
                            middleSp,
                            cache.topDoublets.subset(cache.sortedTops));
  } else {
    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.range(), middleSp,
                            cache.topDoublets.range());
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

  // the cursors live and die with this call, so every middle space point
  // starts its sweep at the front of the subset again
  std::array<SpacePointContainer::ConstSubset, 1> bottomSpGroups{bottomSps};
  std::array<SpacePointContainer::ConstSubset, 1> topSpGroups{topSps};

  createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                            tripletFinder, filter, spacePoints, bottomSpGroups,
                            middleSp, topSpGroups,
                            bottomFinder.config().spacePointsSortedByRadius,
                            bottomFinder.config().deltaRMax,
                            topFinder.config().deltaRMin, outputSeeds);
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
  const float bottomDeltaRMax = bottomFinder.config().deltaRMax;
  const float topDeltaRMin = topFinder.config().deltaRMin;

  if (middleSpGroup.empty()) {
    return;
  }

  // the cursors are ours, so the caller's groups stay untouched and can be
  // handed to another call unchanged
  cache.bottomSpGroups.assign(bottomSpGroups.begin(), bottomSpGroups.end());
  cache.topSpGroups.assign(topSpGroups.begin(), topSpGroups.end());

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

  for (ConstSpacePointProxy spM : middleSpGroup) {
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
    }

    createSeedsFromGroupsImpl(
        *m_logger, cache, bottomFinder, topFinder, tripletFinder, filter,
        spacePoints, cache.bottomSpGroups, spM, cache.topSpGroups,
        spacePointsSortedByRadius, bottomDeltaRMax, topDeltaRMin, outputSeeds);
  }
}

}  // namespace Acts
