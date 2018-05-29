// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
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
#include "Acts/Digitization/Clusterization.hpp"
#include "Acts/Digitization/DigitizationCell.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(merge_clusters1)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    testCells1.push_back(Acts::DigitizationCell(2, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(2, 4, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 6, 1));
    testCells1.push_back(Acts::DigitizationCell(4, 2, 1));
    testCells1.push_back(Acts::DigitizationCell(5, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(9, 5, 1));
    testCells1.push_back(Acts::DigitizationCell(5, 4, 1));
    testCells1.push_back(Acts::DigitizationCell(5, 5, 1));
    testCells1.push_back(Acts::DigitizationCell(6, 5, 1));

    // common edge
    auto mergedCells1 = Acts::createClusters(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells1.size(), 5);
    // common corner
    auto mergedCells2 = Acts::createClusters(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells2.size(), 2);

    std::vector<Acts::DigitizationCell> testCells2;
    testCells2.push_back(Acts::DigitizationCell(2, 3, 1));
    auto mergedCells3 = Acts::createClusters(testCells2, true);

    // test in case just one cell is handed over
    std::vector<Acts::DigitizationCell> testCells3;
    testCells3.push_back(Acts::DigitizationCell(2, 3, 1));
    auto mergedCells4 = Acts::createClusters(testCells3, false);
    BOOST_CHECK_EQUAL(mergedCells4.size(), 1);
    // test in case no cell is handed over
    std::vector<Acts::DigitizationCell> testCells4;
    auto mergedCells5 = Acts::createClusters(testCells4, false);
    BOOST_CHECK_EQUAL(mergedCells5.size(), 0);
  }

  BOOST_AUTO_TEST_CASE(merge_clusters2)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    size_t                              nCells = 99;
    size_t                              delta  = 3;
    size_t nClustersNoTouch = (nCells / delta) * (nCells / delta);

    // Create clusters which are separated by one cell always
    // Clusters have the following shape:
    // -----
    // --##-
    // ---#-
    // -----
    for (size_t i = 0; i < nCells; i += delta) {
      for (size_t j = 0; j < nCells; j += delta) {
        testCells1.push_back(Acts::DigitizationCell(i, j, 1));
        testCells1.push_back(Acts::DigitizationCell(i, j + 1, 1));
        testCells1.push_back(Acts::DigitizationCell(i + 1, j + 1, 1));
      }
    }

    // in a 99x99 grid we get 33x33=1089 clusters
    // common edge
    auto mergedCells1 = Acts::createClusters(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells1.size(), nClustersNoTouch);
    // common corner
    auto mergedCells2 = Acts::createClusters(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells2.size(), nClustersNoTouch);

    // now we add some cells which lead to merging only for common corner (and
    // create new clusters for edge case)
    size_t delta2           = 9;
    size_t nCornerCells     = nCells / delta2 * nCells / delta2;
    size_t nClusters_merged = nClustersNoTouch - nCornerCells * 2;

    for (size_t i = 2; i < nCells; i += delta2) {
      for (size_t j = 2; j < nCells; j += delta2) {
        testCells1.push_back(Acts::DigitizationCell(i, j, 1));
      }
    }

    // common edge
    auto mergedCells3 = Acts::createClusters(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells3.size(), nClustersNoTouch + nCornerCells);

    // common corner
    auto mergedCells4 = Acts::createClusters(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells4.size(), nClusters_merged);

    // now we add some cells which lead to merging also for edge case
    for (size_t i = 2; i < nCells; i += delta2) {
      for (size_t j = 2; j < nCells; j += delta2) {
        testCells1.push_back(Acts::DigitizationCell(i, j - 1, 1));
        testCells1.push_back(Acts::DigitizationCell(i + 1, j, 1));
      }
    }

    // common edge
    auto mergedCells5 = Acts::createClusters(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells5.size(), nClusters_merged);

    // common corner
    auto mergedCells6 = Acts::createClusters(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells6.size(), nClusters_merged);
  }

  BOOST_AUTO_TEST_CASE(merge_cells1)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    size_t                              nCells = 99;
    size_t                              delta  = 3;
    size_t nCellsInClusters = (nCells / delta) * (nCells / delta) * 3;
    // Create clusters which are separated by one cell always
    // Clusters have the following shape:
    // -----
    // --##-
    // ---#-
    // -----
    for (size_t i = 0; i < nCells; i += delta) {
      for (size_t j = 0; j < nCells; j += delta) {
        testCells1.push_back(Acts::DigitizationCell(i, j, 1));
        testCells1.push_back(Acts::DigitizationCell(i, j + 1, 1));
        testCells1.push_back(Acts::DigitizationCell(i + 1, j + 1, 1));
      }
    }

    // digital readout
    auto mergedCells1 = Acts::mergeCells(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells1.size(), nCellsInClusters);
    float data1 = 0;
    for (auto& i : mergedCells1) {
      data1 += i.data;
    }
    BOOST_CHECK_EQUAL(data1, nCellsInClusters);
    // analgogue readout
    auto mergedCells2 = Acts::mergeCells(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells2.size(), nCellsInClusters);
    float data2 = 0;
    for (auto& i : mergedCells2) {
      data2 += i.data;
    }
    BOOST_CHECK_EQUAL(data2, nCellsInClusters);

    // now adding the same cells again - # cells should stay the same after
    // merging
    for (size_t i = 0; i < nCells; i += delta) {
      for (size_t j = 0; j < nCells; j += delta) {
        testCells1.push_back(Acts::DigitizationCell(i, j, 1));
        testCells1.push_back(Acts::DigitizationCell(i, j + 1, 1));
        testCells1.push_back(Acts::DigitizationCell(i + 1, j + 1, 1));
      }
    }

    // digital readout
    auto mergedCells3 = Acts::mergeCells(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells3.size(), nCellsInClusters);
    float data3 = 0;
    for (auto& i : mergedCells3) {
      data3 += i.data;
    }
    BOOST_CHECK_EQUAL(data3, nCellsInClusters);

    // analogue readout - we should have double the energy now
    auto mergedCells4 = Acts::mergeCells(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells4.size(), nCellsInClusters);
    float data4 = 0;
    for (auto& i : mergedCells4) {
      data4 += i.data;
    }
    BOOST_CHECK_EQUAL(data4, nCellsInClusters * 2);
  }

  BOOST_AUTO_TEST_CASE(merge_cells2)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    testCells1.push_back(Acts::DigitizationCell(2, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(2, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 9, 1));
    testCells1.push_back(Acts::DigitizationCell(5, 5, 1));
    testCells1.push_back(Acts::DigitizationCell(5, 6, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 9, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(4, 2, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 9, 1));

    // digital
    auto mergedCells1 = Acts::mergeCells(testCells1, false);
    BOOST_CHECK_EQUAL(mergedCells1.size(), 6);
    float data1 = 0;
    for (auto& i : mergedCells1) {
      data1 += i.data;
    }
    BOOST_CHECK_EQUAL(data1, 6);
    // analogue
    auto mergedCells2 = Acts::mergeCells(testCells1, true);
    BOOST_CHECK_EQUAL(mergedCells2.size(), 6);
    float data2 = 0;
    for (auto& i : mergedCells2) {
      data2 += i.data;
    }
    BOOST_CHECK_EQUAL(data2, 12);

    // test in case only one cell is handed over
    std::vector<Acts::DigitizationCell> testCells2;
    testCells2.push_back(Acts::DigitizationCell(2, 3, 1));
    auto mergedCells3 = Acts::mergeCells(testCells2, true);
    BOOST_CHECK_EQUAL(mergedCells3.size(), 1);

    // test in case no cell is handed over
    std::vector<Acts::DigitizationCell> testCells3;
    auto mergedCells4 = Acts::mergeCells(testCells3, true);
    BOOST_CHECK_EQUAL(mergedCells4.size(), 0);
  }
}
}
