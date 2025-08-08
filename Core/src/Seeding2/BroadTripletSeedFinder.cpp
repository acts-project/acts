// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/BroadTripletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"

#include <Eigen/Dense>

namespace Acts::Experimental {

namespace {

template <typename DoubletCollections>
void createAndFilterTriplets(float rMaxSeedConf,
                             BroadTripletSeedFinder::State& state,
                             BroadTripletSeedFinder::Cache& cache,
                             const TripletSeedFinder& tripletFinder,
                             const BroadTripletSeedFilter& filter,
                             const SpacePointContainer2& spacePoints,
                             DoubletCollections bottomDoublets,
                             const ConstSpacePointProxy2& spM,
                             DoubletCollections topDoublets) {
  for (auto bottomDoublet : bottomDoublets) {
    if (topDoublets.empty()) {
      break;
    }

    cache.tripletTopCandidates.clear();
    tripletFinder.createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                             topDoublets,
                                             cache.tripletTopCandidates);

    auto spB = spacePoints[bottomDoublet.spacePointIndex()];
    // minimum number of compatible top SPs to trigger the filter for a certain
    // middle bottom pair if seedConfirmation is false we always ask for at
    // least one compatible top to trigger the filter
    std::size_t minCompatibleTopSPs = 2;
    if (!filter.config().seedConfirmation || spB.r() > rMaxSeedConf) {
      minCompatibleTopSPs = 1;
    }
    if (filter.config().seedConfirmation &&
        cache.candidatesCollector.nHighQualityCandidates() > 0) {
      minCompatibleTopSPs++;
    }
    // continue if number of top SPs is smaller than minimum required for filter
    if (cache.tripletTopCandidates.size() < minCompatibleTopSPs) {
      continue;
    }

    float zOrigin = spM.z() - spM.r() * bottomDoublet.cotTheta();
    filter.filter2SpFixed(state.filter, cache.filter, spacePoints,
                          bottomDoublet.spacePointIndex(), spM.index(),
                          cache.tripletTopCandidates.topSpacePoints(),
                          cache.tripletTopCandidates.curvatures(),
                          cache.tripletTopCandidates.impactParameters(),
                          zOrigin, cache.candidatesCollector);
  }
}

template <typename SpacePointCollections>
void createSeedsFromGroupsImpl(
    const Logger& logger, BroadTripletSeedFinder::State& state,
    BroadTripletSeedFinder::Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const BroadTripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    SpacePointCollections& bottomSpGroups,
    const ConstSpacePointProxy2& middleSp, SpacePointCollections& topSpGroups,
    SeedContainer2& outputSeeds) {
  MiddleSpInfo middleSpInfo = DoubletSeedFinder::computeMiddleSpInfo(middleSp);

  // create middle-top doublets
  cache.topDoublets.clear();
  for (auto& topSpGroup : topSpGroups) {
    topFinder.createDoublets(middleSp, middleSpInfo, topSpGroup,
                             cache.topDoublets);
  }

  // no top SP found -> cannot form any triplet
  if (cache.topDoublets.empty()) {
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  // apply cut on the number of top SP if seedConfirmation is true
  float rMaxSeedConf = 0;
  if (filter.config().seedConfirmation) {
    // check if middle SP is in the central or forward region
    const bool isForwardRegion =
        middleSp.z() >
            filter.config().centralSeedConfirmationRange.zMaxSeedConf ||
        middleSp.z() <
            filter.config().centralSeedConfirmationRange.zMinSeedConf;
    SeedConfirmationRangeConfig seedConfRange =
        isForwardRegion ? filter.config().forwardSeedConfirmationRange
                        : filter.config().centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the middle SP is
    // in the central or forward region
    std::size_t nTopSeedConf = middleSp.r() > seedConfRange.rMaxSeedConf
                                   ? seedConfRange.nTopForLargeR
                                   : seedConfRange.nTopForSmallR;
    // set max bottom radius for seed confirmation
    rMaxSeedConf = seedConfRange.rMaxSeedConf;
    // continue if number of top SPs is smaller than minimum
    if (cache.topDoublets.size() < nTopSeedConf) {
      ACTS_VERBOSE("Number of top SPs is "
                   << cache.topDoublets.size()
                   << " and is smaller than minimum, returning");
      return;
    }
  }

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  for (auto& bottomSpGroup : bottomSpGroups) {
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

    createAndFilterTriplets(
        rMaxSeedConf, state, cache, tripletFinder, filter, spacePoints,
        cache.bottomDoublets.subset(cache.sortedBottoms), middleSp,
        cache.topDoublets.subset(cache.sortedTops));
  } else {
    createAndFilterTriplets(rMaxSeedConf, state, cache, tripletFinder, filter,
                            spacePoints, cache.bottomDoublets.range(), middleSp,
                            cache.topDoublets.range());
  }

  // retrieve all candidates
  // this collection is already sorted, higher weights first
  const std::size_t numQualitySeeds =
      cache.candidatesCollector.nHighQualityCandidates();
  cache.candidatesCollector.toSortedCandidates(spacePoints,
                                               cache.sortedCandidates);
  filter.filter1SpFixed(state.filter, spacePoints, cache.sortedCandidates,
                        numQualitySeeds, outputSeeds);
}

}  // namespace

BroadTripletSeedFinder::BroadTripletSeedFinder(
    std::unique_ptr<const Logger> logger_)
    : m_logger(std::move(logger_)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument(
        "BroadTripletSeedFinder: logger cannot be null");
  }
}

void BroadTripletSeedFinder::createSeedsFromGroup(
    State& state, Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const BroadTripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    SpacePointContainer2::ConstSubset& bottomSps,
    const ConstSpacePointProxy2& middleSp,
    SpacePointContainer2::ConstSubset& topSps,
    SeedContainer2& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");

  cache.candidatesCollector.setMaxElements(
      filter.config().maxSeedsPerSpMConf,
      filter.config().maxQualitySeedsPerSpMConf);

  std::array<SpacePointContainer2::ConstSubset, 1> bottomSpGroups{bottomSps};
  std::array<SpacePointContainer2::ConstSubset, 1> topSpGroups{topSps};

  createSeedsFromGroupsImpl(*m_logger, state, cache, bottomFinder, topFinder,
                            tripletFinder, filter, spacePoints, bottomSpGroups,
                            middleSp, topSpGroups, outputSeeds);
}

void BroadTripletSeedFinder::createSeedsFromGroups(
    State& state, Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const BroadTripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    const std::span<SpacePointContainer2::ConstRange>& bottomSpGroups,
    const SpacePointContainer2::ConstRange& middleSpGroup,
    const std::span<SpacePointContainer2::ConstRange>& topSpGroups,
    const std::pair<float, float>& radiusRangeForMiddle,
    SeedContainer2& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");
  const bool spacePointsSortedByRadius =
      bottomFinder.config().spacePointsSortedByRadius;

  if (middleSpGroup.empty()) {
    return;
  }

  // initialize cache
  cache.candidatesCollector.setMaxElements(
      filter.config().maxSeedsPerSpMConf,
      filter.config().maxQualitySeedsPerSpMConf);

  if (spacePointsSortedByRadius) {
    // Initialize initial offsets for bottom and top space points with binary
    // search. This requires at least one middle space point to be present which
    // is already checked above.
    const ConstSpacePointProxy2 firstMiddleSp = middleSpGroup.front();
    const float firstMiddleSpR = firstMiddleSp.r();

    for (auto& bottomSpGroup : bottomSpGroups) {
      // Find the first bottom space point that is within the deltaRMax of the
      // first middle space point.
      auto low = std::ranges::lower_bound(
          bottomSpGroup, firstMiddleSpR - bottomFinder.config().deltaRMax, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.r(); });
      bottomSpGroup = bottomSpGroup.subrange(low - bottomSpGroup.begin());
    }

    for (auto& topSpGroup : topSpGroups) {
      // Find the first top space point that is within the deltaRMin of the
      // first middle space point.
      auto low = std::ranges::lower_bound(
          topSpGroup, firstMiddleSpR + topFinder.config().deltaRMin, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.r(); });
      topSpGroup = topSpGroup.subrange(low - topSpGroup.begin());
    }
  }

  for (ConstSpacePointProxy2 spM : middleSpGroup) {
    const float rM = spM.r();

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

    createSeedsFromGroupsImpl(*m_logger, state, cache, bottomFinder, topFinder,
                              tripletFinder, filter, spacePoints,
                              bottomSpGroups, spM, topSpGroups, outputSeeds);
  }
}

}  // namespace Acts::Experimental
