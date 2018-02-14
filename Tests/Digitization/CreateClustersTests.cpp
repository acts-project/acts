// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE CreateClusters Tests

#include <boost/test/included/unit_test.hpp>
// leave blank as
#include <algorithm>
#include <boost/test/data/test_case.hpp>
#include <map>
#include <utility>
#include <vector>
#include "ACTS/Digitization/Digitization.hpp"
#include "ACTS/Digitization/DigitizationCell.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(merge_clusters)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    testCells1.push_back(Acts::DigitizationCell(2, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(2, 4, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 9, 1));

    // common edge
    auto mergedCells1 = Acts::createClusters(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells1.size(), 3);
    // common corner
    auto mergedCells2 = Acts::createClusters(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells1.size(), 2);
  }
}
}
