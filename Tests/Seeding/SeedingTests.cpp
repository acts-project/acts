// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#define BOOST_TEST_MODULE Seeding Tests
#include <boost/test/included/unit_test.hpp>

#include "ACTS/Seeding/BarrelHelixSeedFinder.hpp"
#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Seeding/TrackSeed.hpp"
#include "ACTS/Seeding/detail/SubRange.hpp"

using std::cout;
using namespace Acts::Seeding;

// construct barrel layer w/ n points equidistant in phi
BarrelSpacePoints<size_t>
makeBarrel(double radius, int nPoints)
{
  BarrelSpacePoints<size_t> barrel(radius);

  for (int i = 0; i < nPoints; ++i) {
    double phi = 2 * M_PI * i / nPoints;
    barrel.points.emplace_back(
        radius * std::cos(phi), radius * std::sin(phi), 0, i);
  }
  barrel.sortByPhi();
  return barrel;
}

template <typename Identifier>
void
print(std::ostream& os, const BarrelSpacePoints<Identifier>& barrel)
{
  for (const auto& point : barrel.points) os << point << '\n';
}

template <typename Identifier, size_t N>
void
print(std::ostream& os, const std::vector<TrackSeed<Identifier, N>>& seeds)
{
  for (const auto& seed : seeds) {
    os << seed << '\n';
    for (size_t i = 0; i < seed.kNumPoints; ++i) {
      os << "  " << seed.point(i) << '\n';
    }
  }
}

BOOST_AUTO_TEST_CASE(PhiRangeTest)
{
  using detail::makeRangePhi;

  auto layer = makeBarrel(10, 16);
  auto dphi  = 2 * M_PI / 16;

  print(std::cout, layer);

  {
    // empty range, linear
    auto range = layer.rangePhiDelta(0.5 * dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(range.size(), 0);
    BOOST_CHECK(range.begin() == range.end());
  }
  {
    // empty range at edge
    auto range = layer.rangePhiDelta(2 * M_PI - 0.5 * dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(range.size(), 0);
    BOOST_CHECK(range.begin() == range.end());
  }
  {
    // linear range, no wrap-around
    auto range = layer.rangePhiDelta(1 * dphi, 1 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(range.size(), 3);
    BOOST_CHECK_CLOSE(it->phi(), 0 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), 1 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), 2 * dphi, 1e-6);
    ++it;
    // BOOST_CHECK_CLOSE(it->phi(), 3 * dphi, 1e-6);
    // ++it;
    BOOST_CHECK(it == range.end());
  }
  {
    // wrap around at edge
    auto range = layer.rangePhiDelta(M_PI, 2 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(range.size(), 5);
    BOOST_CHECK_CLOSE(it->phi(), M_PI - 2 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), M_PI - 1 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), M_PI - 0 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), -M_PI + 1 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE(it->phi(), -M_PI + 2 * dphi, 1e-6);
    ++it;
    BOOST_CHECK(it == range.end());
  }
}

BOOST_AUTO_TEST_CASE(BarrelHelixSeedFinderTest)
{
  size_t pointsPerLayer = 16;
  double dphi           = 2 * M_PI / pointsPerLayer;

  auto layer0 = makeBarrel(10, pointsPerLayer);
  auto layer1 = makeBarrel(30, pointsPerLayer);
  auto layer2 = makeBarrel(50, pointsPerLayer);

  print(std::cout, layer0);
  print(std::cout, layer1);
  print(std::cout, layer2);

  // hard cuts, only straight tracks
  {
    HelixSeedConfig cfg;
    cfg.deltaPhi01 = 0.1 * dphi;
    cfg.deltaPhi12 = 0.1 * dphi;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);
    print(std::cout, seeds);

    BOOST_CHECK_EQUAL(seeds.size(), pointsPerLayer);
  }
  // looser cuts, multiple combinations for each inner hit
  {
    HelixSeedConfig cfg;
    cfg.deltaPhi01 = 1.1 * dphi;
    cfg.deltaPhi12 = 1.1 * dphi;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);
    print(std::cout, seeds);

    BOOST_CHECK_GT(seeds.size(), pointsPerLayer);
  }
}
