// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#define BOOST_TEST_MODULE Seeding Tools Tests
#include <boost/test/included/unit_test.hpp>

#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Seeding/detail/geometry.hpp"
#include "ACTS/Seeding/detail/cyclic_range.hpp"
#include "SeedingTestsCommon.hpp"

BOOST_AUTO_TEST_CASE(PhiRangeTest)
{
  using Acts::detail::makeRangePhi;

  auto layer              = makeBarrel(10, 16);
  auto dphi               = 2 * M_PI / 16;
  auto compareSpacePoints = [](const auto& a, const auto& b) {
    return (a.identifier() == b.identifier());
  };

  print(layer);

  {
    // empty range, in the middle
    auto range = layer.rangeDeltaPhi(dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  {
    // empty range at left edge
    auto range = layer.rangeDeltaPhi(0, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  {
    // empty range at right edge
    auto range = layer.rangeDeltaPhi(M_PI, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  {
    // full linear range
    auto range = layer.rangeDeltaPhi(0, M_PI - 0.1 * dphi);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()),
                      layer.points.size());
    BOOST_CHECK(range.begin() != range.end());
    BOOST_CHECK(std::equal(
        range.begin(), range.end(), layer.points.begin(), compareSpacePoints));
  }
  {
    // linear range, no wrap-around
    auto range = layer.rangeDeltaPhi(1.5 * dphi, 1.1 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 3);
    BOOST_CHECK_CLOSE((*it).phi(), 0.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), 1.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), 2.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK(!(it != range.end()));
  }
  {
    // wrap around at edge
    auto range = layer.rangeDeltaPhi(M_PI, 3.1 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 6);
    BOOST_CHECK_CLOSE((*it).phi(), M_PI - 2.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), M_PI - 1.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), M_PI - 0.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), -M_PI + 0.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), -M_PI + 1.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK_CLOSE((*it).phi(), -M_PI + 2.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK(!(it != range.end()));
  }
}

BOOST_AUTO_TEST_CASE(ZeroCurvatureCircleTest)
{
  using Acts::detail::calcCircleCurvature;

  Acts::Vector3D delta(1, 1, 0);

  BOOST_CHECK_CLOSE(calcCircleCurvature(delta, delta), 0, 1e-6);
  BOOST_CHECK_CLOSE(calcCircleCurvature(-delta, -delta), 0, 1e-6);
}

BOOST_AUTO_TEST_CASE(SignCurvatureCircleTest)
{
  using Acts::detail::calcCircleCurvature;

  Acts::Vector3D d01(1, 0, 0);
  Acts::Vector3D d12pos(1, 1, 0);
  Acts::Vector3D d12neg(1, -1, 0);

  BOOST_CHECK_LT(calcCircleCurvature(d01, d12pos), 0);
  BOOST_CHECK_GT(calcCircleCurvature(d01, d12neg), 0);
}

BOOST_AUTO_TEST_CASE(EstimateD0CircleTest)
{
  using Acts::Vector3D;
  using Acts::detail::estimateCircleD0;

  // track position on the y-axis pointing towards positive x
  double phi = 0;
  // vanishing curvature
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, 1, 0), phi, 0), 1, 1e-6);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(10, 1, 0), phi, 0), 1, 1e-6);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, -1, 0), phi, 0), -1, 1e-6);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(10, -1, 0), phi, 0), -1, 1e-6);
  // small curvature (epsilon is given in percent, *not* fraction)
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, 1, 0), phi, 1e-4), 1, 0.1);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, 1, 0), phi, -1e-4), 1, 0.1);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, -1, 0), phi, 1e-4), -1, 0.1);
  BOOST_CHECK_CLOSE(estimateCircleD0(Vector3D(0, -1, 0), phi, -1e-4), -1, 0.1);
}
