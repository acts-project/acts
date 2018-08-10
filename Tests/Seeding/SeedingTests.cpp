// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#define BOOST_TEST_MODULE Seeding Tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/Seeding/BarrelSeedFinder.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/TrackSeed.hpp"
#include "SeedingTestsCommon.hpp"

using namespace Acts::Seeding;

BOOST_AUTO_TEST_CASE(HelixBarrelSeedFinderTest)
{
  size_t pointsPerLayer = 16;
  double dphi           = 2 * M_PI / pointsPerLayer;

  auto layer0 = makeBarrel(10, pointsPerLayer);
  auto layer1 = makeBarrel(30, pointsPerLayer);
  auto layer2 = makeBarrel(50, pointsPerLayer);

  // hard cuts, only straight tracks
  {
    HelixSeedConfig cfg;
    cfg.rangePhi1     = 0.001 * dphi;
    cfg.rangePhi2     = 0.001 * dphi;
    cfg.maxDeltaTheta = 0.2;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);

    BOOST_CHECK_EQUAL(seeds.size(), pointsPerLayer);
  }
  // medium cuts, hit1 is matched, hit2 pickups 3 hits
  {
    HelixSeedConfig cfg;
    cfg.rangePhi1     = 0.001 * dphi;
    cfg.rangePhi2     = 1.001 * dphi;  // matching hit +- 1 neighbor
    cfg.maxDeltaTheta = 0.2;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);

    BOOST_CHECK_EQUAL(seeds.size(), 3 * pointsPerLayer);
  }
  // loose cuts, hit1 picks up 3 hits, each hit2 picks up 3 hits
  {
    HelixSeedConfig cfg;
    cfg.rangePhi1     = 1.001 * dphi;  // matching hit +- 1 neighbor
    cfg.rangePhi2     = 1.001 * dphi;  // matching hit +- 1 neighbor
    cfg.maxDeltaTheta = 0.2;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);

    BOOST_CHECK_EQUAL(seeds.size(), 9 * pointsPerLayer);
  }
}
