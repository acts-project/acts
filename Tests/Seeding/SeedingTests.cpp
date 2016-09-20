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

#include "ACTS/Seeding/BarrelSeedFinder.hpp"
#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Seeding/TrackSeed.hpp"

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
  barrel.sort();
  return barrel;
}

template <typename Identifier>
void
print(const BarrelSpacePoints<Identifier>& barrel)
{
  for (const auto& point : barrel.points) std::cout << point << '\n';
}

template <typename Identifier, size_t N>
void
print(const std::vector<TrackSeed<Identifier, N>>& seeds)
{
  for (const auto& seed : seeds) seed.print(std::cout);
}

BOOST_AUTO_TEST_CASE(PhiRangeTest)
{
  using detail::makeRangePhi;

  auto layer = makeBarrel(10, 16);
  auto dphi  = 2 * M_PI / 16;

  print(layer);

  {
    // empty range, linear
    auto range = layer.rangeDeltaPhi(0.5 * dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(range.size(), 0);
    BOOST_CHECK(range.begin() == range.end());
  }
  {
    // empty range at edge
    auto range = layer.rangeDeltaPhi(2 * M_PI - 0.5 * dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(range.size(), 0);
    BOOST_CHECK(range.begin() == range.end());
  }
  {
    // linear range, no wrap-around
    auto range = layer.rangeDeltaPhi(1 * dphi, 1 * dphi);
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
    auto range = layer.rangeDeltaPhi(M_PI, 2 * dphi);
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

BOOST_AUTO_TEST_CASE(SignCircleCurvatureTest)
{
  using Acts::Vector3D;

  Vector3D p0(0, 0, 0);
  Vector3D p1(1, 0, 0);
  Vector3D p2pos(2, 1, 0);
  Vector3D p2neg(2, -1, 0);

  BOOST_CHECK_LT(detail::calcCircleCurvature(p0, p1, p2pos), 0);
  BOOST_CHECK_GT(detail::calcCircleCurvature(p0, p1, p2neg), 0);
}

BOOST_AUTO_TEST_CASE(HelixBarrelSeedFinderTest)
{
  size_t pointsPerLayer = 16;
  double dphi           = 2 * M_PI / pointsPerLayer;

  auto layer0 = makeBarrel(10, pointsPerLayer);
  auto layer1 = makeBarrel(30, pointsPerLayer);
  auto layer2 = makeBarrel(50, pointsPerLayer);

  print(layer0);
  print(layer1);
  print(layer2);

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
    cfg.rangePhi2     = 1.2 * dphi;    // matching hit +- 1 neighbor
    cfg.maxDeltaTheta = 0.2;
    TrackSeeds3<size_t> seeds;

    findHelixSeeds(cfg, layer0, layer1, layer2, seeds);
    print(seeds);

    BOOST_CHECK_EQUAL(seeds.size(), 9 * pointsPerLayer);
  }
}
