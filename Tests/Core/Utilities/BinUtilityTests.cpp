// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE BinUtility Tests
#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"

namespace Acts {
namespace Test {

  namespace tt = boost::test_tools;

  // OPEN - equidistant binning tests
  BOOST_AUTO_TEST_CASE(BinUtility_equidistant_binning)
  {
    Vector3D xyzPosition(1.5, 2.5, 3.5);
    Vector3D edgePosition(0.5, 0.5, 0.5);

    // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
    BinUtility xUtil_eq(10, 0., 10., open, binX);
    BinUtility yUtil_eq(10, 0., 10., open, binY);
    BinUtility zUtil_eq(10, 0., 10., open, binZ);
    BOOST_CHECK_EQUAL(xUtil_eq.bins(), (size_t)10);
    // make it 2-dim
    BinUtility xyUtil_eq(10, 0., 10., open, binX);
    xyUtil_eq += yUtil_eq;
    BOOST_CHECK_EQUAL(xyUtil_eq.bins(), 100);
    // make it 3-dim
    BinUtility xyzUtil_eq(xyUtil_eq);
    xyzUtil_eq += zUtil_eq;
    BOOST_CHECK_EQUAL(xyzUtil_eq.bins(), 1000);
    // check the dimensions
    BOOST_CHECK_EQUAL(xUtil_eq.dimensions(), 1);
    BOOST_CHECK_EQUAL(xyUtil_eq.dimensions(), 2);
    BOOST_CHECK_EQUAL(xyzUtil_eq.dimensions(), 3);

    // bin triples and clusters
    auto xTriple   = xUtil_eq.binTriple(xyzPosition);
    auto xyTriple  = xyUtil_eq.binTriple(xyzPosition);
    auto xyzTriple = xyzUtil_eq.binTriple(xyzPosition);

    BOOST_CHECK_EQUAL(xTriple[0], 1);
    BOOST_CHECK_EQUAL(xTriple[1], 0);
    BOOST_CHECK_EQUAL(xTriple[2], 0);

    BOOST_CHECK_EQUAL(xyTriple[0], 1);
    BOOST_CHECK_EQUAL(xyTriple[1], 2);
    BOOST_CHECK_EQUAL(xyTriple[2], 0);

    BOOST_CHECK_EQUAL(xyzTriple[0], 1);
    BOOST_CHECK_EQUAL(xyzTriple[1], 2);
    BOOST_CHECK_EQUAL(xyzTriple[2], 3);

    // Full range
    std::vector<size_t> xRangeCheck0   = {0, 1, 2};
    std::vector<size_t> xyRangeCheck1  = {1, 2, 3};
    std::vector<size_t> xyzRangeCheck2 = {2, 3, 4};

    auto xNrange0 = xUtil_eq.neighbourRange(xyzPosition, 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(xNrange0.begin(),
                                  xNrange0.end(),
                                  xRangeCheck0.begin(),
                                  xRangeCheck0.end());

    auto xNrange1 = xUtil_eq.neighbourRange(xyzPosition, 1);
    BOOST_CHECK_EQUAL(xNrange1.size(), 1);
    BOOST_CHECK_EQUAL(xNrange1[0], 0);

    auto xNrange2 = xUtil_eq.neighbourRange(xyzPosition, 2);
    BOOST_CHECK_EQUAL(xNrange2.size(), 1);
    BOOST_CHECK_EQUAL(xNrange2[0], 0);

    auto xyNrange1 = xyUtil_eq.neighbourRange(xyzPosition, 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(xyNrange1.begin(),
                                  xyNrange1.end(),
                                  xyRangeCheck1.begin(),
                                  xyRangeCheck1.end());

    auto xyzNrange2 = xyzUtil_eq.neighbourRange(xyzPosition, 2);
    BOOST_CHECK_EQUAL_COLLECTIONS(xyzNrange2.begin(),
                                  xyzNrange2.end(),
                                  xyzRangeCheck2.begin(),
                                  xyzRangeCheck2.end());

    // Partial range
    std::vector<size_t> xEdgeCheck = {0, 1};
    auto                xEdgeRange = xUtil_eq.neighbourRange(edgePosition, 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(xEdgeRange.begin(),
                                  xEdgeRange.end(),
                                  xEdgeCheck.begin(),
                                  xEdgeCheck.end());
  }
}
}
