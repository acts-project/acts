// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <limits>

using namespace Acts;

namespace ActsTests {

namespace {

// filterTripletsMiddleFixed reads spacePoints[index].index() for every
// candidate pushed into the collector, so the container just needs enough
// entries for the indices used below; none of the coordinate columns are
// read by the filter.
SpacePointContainer makeSpacePoints(std::size_t n) {
  SpacePointContainer spacePoints;
  spacePoints.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    spacePoints.createSpacePoint();
  }
  return spacePoints;
}

// Pushes `n` low-quality candidates directly into the filter's collector,
// bypassing seed-confirmation push logic (filterTripletTopCandidates),
// which is not exercised by this test suite. Weights are strictly
// descending and every candidate uses distinct bottom/top space points, so
// no candidate can be skipped by the best-seed-quality bookkeeping in
// filterTripletsMiddleFixed (see getBestSeedQuality/setBestSeedQuality in
// BroadTripletSeedFilter.cpp): with `nLow` the collector's low-quality
// capacity, the number of seeds filterTripletsMiddleFixed emits is exactly
// min(n, nLow, maxSeedsPerSpM + 1).
void pushLowQualityCandidates(BroadTripletSeedFilter::State& state,
                              std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    const auto index = static_cast<SpacePointIndex>(i);
    const float weight = static_cast<float>(n - i);
    state.candidatesCollector.push(/*spB=*/index, /*spM=*/0, /*spT=*/index,
                                   weight, /*zOrigin=*/0.f,
                                   /*isQuality=*/false);
  }
}

std::size_t runFilter(const BroadTripletSeedFilter::Config& config,
                      std::size_t nCandidates, std::size_t nSpacePoints) {
  BroadTripletSeedFilter::State state;
  BroadTripletSeedFilter::Cache cache;
  auto logger =
      getDefaultLogger("BroadTripletSeedFilterTest", Logging::Level::WARNING);
  BroadTripletSeedFilter filter(config, state, cache, *logger);

  pushLowQualityCandidates(state, nCandidates);

  auto spacePoints = makeSpacePoints(nSpacePoints);
  SeedContainer seeds;
  filter.filterTripletsMiddleFixed(spacePoints, seeds);
  return seeds.size();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(BroadTripletSeedFilterTests)

// Regression test for the collector-sizing bug: with seedConfirmation off,
// the low-quality collector used to be sized from maxSeedsPerSpMConf
// regardless of maxSeedsPerSpM, so a maxSeedsPerSpM above maxSeedsPerSpMConf
// - 1 had no effect. Fails before the fix (yields 5, not 11).
BOOST_AUTO_TEST_CASE(MaxSeedsPerSpMIsNotClampedByConfWhenConfirmationOff) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = 10;
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, /*nCandidates=*/20, /*nSpacePoints=*/20),
                    11u);
}

// The physmon-relevant case: maxSeedsPerSpM below the default
// maxSeedsPerSpMConf already bound before the fix, so this stays unchanged
// (min(20, 5, 2) == min(20, 2, 2) == 2).
BOOST_AUTO_TEST_CASE(MaxSeedsPerSpMBelowDefaultConfIsUnchanged) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = 1;
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, /*nCandidates=*/20, /*nSpacePoints=*/20),
                    2u);
}

// maxSeedsPerSpM is unsigned int, and the collector capacity is a
// std::uint32_t saturating at CandidatesForMiddleSp::kNoSize. Regression
// test for the +1 overflowing to 0 (which would reject every candidate)
// when maxSeedsPerSpM is already the maximum representable value.
BOOST_AUTO_TEST_CASE(MaxSeedsPerSpMAtLimitDoesNotOverflow) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = std::numeric_limits<unsigned int>::max();
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, /*nCandidates=*/20, /*nSpacePoints=*/20),
                    20u);
}

// With seedConfirmation on, the collector is still sized from the
// confirmation parameters, unaffected by the fix. This only pins the
// constructor's sizing decision: candidates are pushed directly into the
// collector here, bypassing filterTripletTopCandidates, so this is not a
// test of seed-confirmation candidate selection.
BOOST_AUTO_TEST_CASE(SeedConfirmationOnKeepsConfSizing) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = true;
  config.maxSeedsPerSpM = 10;
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, /*nCandidates=*/20, /*nSpacePoints=*/20),
                    5u);
}

// Pins the deliberate choice to size the high-quality collection to 0 when
// seedConfirmation is off: nothing in the filter ever pushes a
// high-quality candidate in that mode (only filterTripletTopCandidates
// does, gated on seedConfirmation), so the budget is unreachable through
// the filter's own code paths. This test pushes one directly to observe
// the constructor's sizing decision; reverting the nHigh choice to
// maxQualitySeedsPerSpMConf only requires updating this one test.
BOOST_AUTO_TEST_CASE(HighQualityCollectionIsUnusedWhenConfirmationOff) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = 10;
  config.maxQualitySeedsPerSpMConf = 5;

  BroadTripletSeedFilter::State state;
  BroadTripletSeedFilter::Cache cache;
  auto logger =
      getDefaultLogger("BroadTripletSeedFilterTest", Logging::Level::WARNING);
  BroadTripletSeedFilter filter(config, state, cache, *logger);

  const bool added = state.candidatesCollector.push(
      /*spB=*/0, /*spM=*/1, /*spT=*/2, /*weight=*/1.f, /*zOrigin=*/0.f,
      /*isQuality=*/true);
  BOOST_CHECK(!added);
  BOOST_CHECK_EQUAL(state.candidatesCollector.nHighQualityCandidates(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
