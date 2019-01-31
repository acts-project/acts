// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Seeding Tools Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <iostream>

#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/detail/cyclic_range.hpp"
#include "Acts/Seeding/detail/geometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include "SeedingTestsCommon.hpp"

namespace data = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(ContainerSizeOnePhiCyclicRangeTest)
{
  using Acts::detail::makeRangePhi;

  auto points = makeBarrel(1.0, 1).points;
  auto eps    = 0.001;

  BOOST_TEST_CONTEXT("full range linear")
  {
    auto range = makeRangePhi(points, -M_PI + eps, M_PI - eps);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 0lu);
    BOOST_CHECK_EQUAL(range.end().index, 1lu);
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 1);
  }
  BOOST_TEST_CONTEXT("linear w/o first element")
  {
    auto range = makeRangePhi(points, eps, M_PI - eps);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
  }
  BOOST_TEST_CONTEXT("wrap-around w/o first element")
  {
    auto range = makeRangePhi(points, M_PI - eps, -eps);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
  }
}

BOOST_DATA_TEST_CASE(VariableContainerSizePhiCyclicRangeTest,
                     data::xrange(2, 9),
                     n)
{
  using Acts::detail::makeRangePhi;

  auto points = makeBarrel(1.0, n).points;
  auto dphi   = 2 * M_PI / n;

  // bounds with epsilons to ensure correct float comparisons
  double low  = -M_PI + 0.001 * dphi;
  double high = M_PI - 0.001 * dphi;
  // large enough to (de-)select one element but small enough to still retain
  // ordering after wrap-around
  double delta = 0.9 * dphi;
  BOOST_TEST_CONTEXT("linear range full")
  {
    auto range = makeRangePhi(points, low, high);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 0lu);
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)n);
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), n);
  }
  BOOST_TEST_CONTEXT("linear range w/o first element")
  {
    auto range = makeRangePhi(points, low + delta, high);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 1lu);
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)n);
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), (n - 1));
  }
  BOOST_TEST_CONTEXT("linear range w/o last element")
  {
    auto range = makeRangePhi(points, low, high - delta);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 0ul);
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)(n - 1));
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), (n - 1));
  }
  BOOST_TEST_CONTEXT("wrap-around full range w/ one element after wrap")
  {
    auto range = makeRangePhi(points, low + delta, high + delta);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 1ul);
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)(1 + n));
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), n);
  }
  BOOST_TEST_CONTEXT("wrap-around full range w/ one element before wrap")
  {
    auto range = makeRangePhi(points, low - delta, high - delta);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, (unsigned int)(n - 1));
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)(n - 1 + n));
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), n);
  }
  BOOST_TEST_CONTEXT("wrap-around w/ only last and first element")
  {
    auto range = makeRangePhi(points, high - delta, low + delta);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, (unsigned int)(n - 1));
    BOOST_CHECK_EQUAL(range.end().index, (unsigned int)(n + 1));
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 2);
  }
  BOOST_TEST_CONTEXT("wrap-around boundaries but only 1 element after wrap")
  {
    auto range = makeRangePhi(points, high, low + delta);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(range.begin().index, 0ul);
    BOOST_CHECK_EQUAL(range.end().index, 1ul);
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 1);
  }
}

BOOST_AUTO_TEST_CASE(ExampleContainerPhiCyclicRangeTest)
{
  using Acts::detail::makeRangePhi;

  auto layer           = makeBarrel(10, 16);
  auto dphi            = 2 * M_PI / 16;
  auto compSpacePoints = [](const auto& a, const auto& b) {
    return (a.identifier() == b.identifier());
  };

  BOOST_TEST_CONTEXT("empty range linear")
  {
    auto range = layer.rangeDeltaPhi(dphi, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  BOOST_TEST_CONTEXT("empty range left edge")
  {
    auto range = layer.rangeDeltaPhi(0, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  BOOST_TEST_CONTEXT("empty range right edge")
  {
    auto range = layer.rangeDeltaPhi(M_PI, 0.4 * dphi);
    BOOST_CHECK(range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 0);
    BOOST_CHECK(!(range.begin() != range.end()));
  }
  BOOST_TEST_CONTEXT("full linear range")
  {
    auto range = layer.rangeDeltaPhi(0, M_PI - 0.1 * dphi);
    BOOST_CHECK(!range.empty());
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()),
                      (int)layer.points.size());
    BOOST_CHECK(std::equal(
        range.begin(), range.end(), layer.points.begin(), compSpacePoints));
  }
  BOOST_TEST_CONTEXT("linear range w/o wrap-around")
  {
    auto range = layer.rangeDeltaPhi(1.5 * dphi, 1.1 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 3);
    CHECK_CLOSE_REL((*it).phi(), 0.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), 1.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), 2.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK(!(it != range.end()));
  }
  BOOST_TEST_CONTEXT("partial range w/ wrap-around")
  {
    auto range = layer.rangeDeltaPhi(M_PI, 2.6 * dphi);
    auto it    = range.begin();
    BOOST_CHECK_EQUAL(std::distance(range.begin(), range.end()), 6);
    CHECK_CLOSE_REL((*it).phi(), M_PI - 2.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), M_PI - 1.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), M_PI - 0.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), -M_PI + 0.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), -M_PI + 1.5 * dphi, 1e-6);
    ++it;
    CHECK_CLOSE_REL((*it).phi(), -M_PI + 2.5 * dphi, 1e-6);
    ++it;
    BOOST_CHECK(!(it != range.end()));
  }
}

BOOST_AUTO_TEST_CASE(ZeroCurvatureCircleTest)
{
  using Acts::detail::calcCircleCurvature;

  Acts::Vector3D delta(1, 1, 0);

  CHECK_SMALL(calcCircleCurvature(delta, delta), 1e-9);
  CHECK_SMALL(calcCircleCurvature(-delta, -delta), 1e-9);
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
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, 1, 0), phi, 0), 1, 1e-6);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(10, 1, 0), phi, 0), 1, 1e-6);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, -1, 0), phi, 0), -1, 1e-6);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(10, -1, 0), phi, 0), -1, 1e-6);
  // small curvature
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, 1, 0), phi, 1e-4), 1, 1e-3);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, 1, 0), phi, -1e-4), 1, 1e-3);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, -1, 0), phi, 1e-4), -1, 1e-3);
  CHECK_CLOSE_REL(estimateCircleD0(Vector3D(0, -1, 0), phi, -1e-4), -1, 1e-3);
}
