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

#include <cstdint>
#include <limits>

using namespace Acts;

namespace ActsTests {

namespace {

// Returns the number of seeds kept from `nCandidates` low-quality candidates
std::size_t runFilter(const BroadTripletSeedFilter::Config& config,
                      std::size_t nCandidates) {
  BroadTripletSeedFilter::State state;
  BroadTripletSeedFilter::Cache cache;
  auto logger =
      getDefaultLogger("BroadTripletSeedFilterTest", Logging::Level::WARNING);
  BroadTripletSeedFilter filter(config, state, cache, *logger);

  SpacePointContainer spacePoints;
  for (std::size_t i = 0; i < nCandidates; ++i) {
    spacePoints.createSpacePoint();
    const auto index = static_cast<SpacePointIndex>(i);
    state.candidatesCollector.push(
        index, 0, index, static_cast<float>(nCandidates - i), 0.f, false);
  }

  SeedContainer seeds;
  filter.filterTripletsMiddleFixed(spacePoints, seeds);
  return seeds.size();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(BroadTripletSeedFilterTests)

// maxSeedsPerSpMConf must not limit maxSeedsPerSpM without seed confirmation
BOOST_AUTO_TEST_CASE(MaxSeedsPerSpMIsNotClampedByConfWhenConfirmationOff) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = 10;
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, 20), 11u);
}

// maxSeedsPerSpM + 1 must not wrap around to a capacity of zero
BOOST_AUTO_TEST_CASE(MaxSeedsPerSpMAtLimitDoesNotOverflow) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = false;
  config.maxSeedsPerSpM = std::numeric_limits<std::uint32_t>::max();

  BOOST_CHECK_EQUAL(runFilter(config, 20), 20u);
}

// With seed confirmation the collector is sized from maxSeedsPerSpMConf
BOOST_AUTO_TEST_CASE(SeedConfirmationOnKeepsConfSizing) {
  BroadTripletSeedFilter::Config config;
  config.seedConfirmation = true;
  config.maxSeedsPerSpM = 10;
  config.maxSeedsPerSpMConf = 5;

  BOOST_CHECK_EQUAL(runFilter(config, 20), 5u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
