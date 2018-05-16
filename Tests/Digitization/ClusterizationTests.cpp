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
#include <chrono>
#include <ctime>
#include <map>
#include <utility>
#include <vector>
#include "Acts/Digitization/Clusterization.hpp"
#include "Acts/Digitization/DigitizationCell.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// This test tests the clusterization of cells which belong to the same
  /// cluster for 8-cell/4-cell merging, digital/analogue readout and with a
  /// possible energy cut applied
  /// The grid with cells should cover all different cases:
  ///
  /// 1 0 0 0 2 0 0 0 0 2
  /// 0 2 0 1 0 0 0 2 0 2
  /// 0 0 0 0 0 0 0 1 1 2
  /// 2 0 0 0 0 1 0 0 0 0
  /// 1 0 0 1 0 2 0 0 0 2
  /// 0 0 2 0 0 1 0 1 0 1
  /// 0 0 0 1 1 0 0 2 0 1
  /// 0 0 0 0 0 0 0 0 2 0
  /// 1 2 2 0 0 0 0 0 0 0
  BOOST_AUTO_TEST_CASE(create_Clusters1)
  {
    size_t              nBins0          = 10;
    size_t              nBins1          = 9;
    std::vector<size_t> clusterSizes    = {2, 2, 2, 3, 6, 6, 7};
    std::vector<size_t> clusterSizesCut = {1, 1, 1, 1, 1, 1, 1, 2, 2, 3};
    std::vector<size_t> clusterSizesEdge
        = {1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 6};
    size_t nClusters     = clusterSizes.size();
    size_t nClustersCut  = clusterSizesCut.size();
    size_t nClustersEdge = clusterSizesEdge.size();

    std::vector<Acts::DigitizationCell> testCells;
    // add cells covering all cases
    testCells.push_back(Acts::DigitizationCell(0, 0, 1));
    testCells.push_back(Acts::DigitizationCell(0, 3, 1));
    testCells.push_back(Acts::DigitizationCell(0, 4, 1));
    testCells.push_back(Acts::DigitizationCell(0, 8, 1));
    testCells.push_back(Acts::DigitizationCell(1, 1, 1));
    testCells.push_back(Acts::DigitizationCell(1, 8, 1));
    testCells.push_back(Acts::DigitizationCell(2, 5, 1));
    testCells.push_back(Acts::DigitizationCell(2, 8, 1));
    testCells.push_back(Acts::DigitizationCell(3, 1, 1));
    testCells.push_back(Acts::DigitizationCell(3, 4, 1));
    testCells.push_back(Acts::DigitizationCell(3, 6, 1));
    testCells.push_back(Acts::DigitizationCell(4, 0, 1));
    testCells.push_back(Acts::DigitizationCell(4, 6, 1));
    testCells.push_back(Acts::DigitizationCell(5, 3, 1));
    testCells.push_back(Acts::DigitizationCell(5, 4, 1));
    testCells.push_back(Acts::DigitizationCell(5, 5, 1));
    testCells.push_back(Acts::DigitizationCell(7, 1, 1));
    testCells.push_back(Acts::DigitizationCell(7, 2, 1));
    testCells.push_back(Acts::DigitizationCell(7, 5, 1));
    testCells.push_back(Acts::DigitizationCell(7, 6, 1));
    testCells.push_back(Acts::DigitizationCell(8, 2, 1));
    testCells.push_back(Acts::DigitizationCell(8, 7, 1));
    testCells.push_back(Acts::DigitizationCell(9, 0, 1));
    testCells.push_back(Acts::DigitizationCell(9, 1, 1));
    testCells.push_back(Acts::DigitizationCell(9, 2, 1));
    testCells.push_back(Acts::DigitizationCell(9, 4, 1));
    testCells.push_back(Acts::DigitizationCell(9, 5, 1));
    testCells.push_back(Acts::DigitizationCell(9, 6, 1));

    size_t nCellsWithoutDuplicates = testCells.size();

    // add duplicates
    testCells.push_back(Acts::DigitizationCell(0, 3, 1));
    testCells.push_back(Acts::DigitizationCell(1, 1, 1));
    testCells.push_back(Acts::DigitizationCell(1, 8, 1));
    testCells.push_back(Acts::DigitizationCell(2, 5, 1));
    testCells.push_back(Acts::DigitizationCell(2, 8, 1));
    testCells.push_back(Acts::DigitizationCell(4, 0, 1));
    testCells.push_back(Acts::DigitizationCell(5, 4, 1));
    testCells.push_back(Acts::DigitizationCell(7, 1, 1));
    testCells.push_back(Acts::DigitizationCell(7, 6, 1));
    testCells.push_back(Acts::DigitizationCell(8, 7, 1));
    testCells.push_back(Acts::DigitizationCell(9, 0, 1));
    testCells.push_back(Acts::DigitizationCell(9, 1, 1));
    testCells.push_back(Acts::DigitizationCell(9, 2, 1));
    testCells.push_back(Acts::DigitizationCell(9, 4, 1));

    size_t nCells = testCells.size();

    // Common Corner, digital,no energy cut
    // createClusters
    auto mergedCells1 = Acts::createClusters<Acts::DigitizationCell>(
        testCells, nBins0, nBins1, true, false, 0.);
    // check number of clusters
    BOOST_CHECK_EQUAL(mergedCells1.size(), nClusters);

    float               data1 = 0;
    std::vector<size_t> clusterSizes1;
    for (size_t i = 0; i < mergedCells1.size(); i++) {
      auto cells = mergedCells1.at(i);
      for (auto& cell : cells) {
        data1 += cell.data;
      }
      // check the cluster sizes
      clusterSizes1.push_back(cells.size());
    }
    // check the cluster sizes
    std::sort(clusterSizes1.begin(), clusterSizes1.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizes.begin(),
                                  clusterSizes.end(),
                                  clusterSizes1.begin(),
                                  clusterSizes1.end());
    // check cells
    BOOST_CHECK_CLOSE(data1, nCellsWithoutDuplicates, 10e-5);

    // Common Edge, digital,no energy cut
    // createClusters
    auto mergedCells2 = Acts::createClusters<Acts::DigitizationCell>(
        testCells, nBins0, nBins1, false, false, 0.);
    // check number of clusters
    BOOST_CHECK_EQUAL(mergedCells2.size(), nClustersEdge);

    float               data2 = 0;
    std::vector<size_t> clusterSizes2;
    for (size_t i = 0; i < mergedCells2.size(); i++) {
      auto cells = mergedCells2.at(i);
      for (auto& cell : cells) {
        data2 += cell.data;
      }
      // check the cluster sizes
      clusterSizes2.push_back(cells.size());
    }
    // check the cluster sizes
    std::sort(clusterSizes2.begin(), clusterSizes2.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizesEdge.begin(),
                                  clusterSizesEdge.end(),
                                  clusterSizes2.begin(),
                                  clusterSizes2.end());
    // check cells
    BOOST_CHECK_CLOSE(data2, nCellsWithoutDuplicates, 10e-5);

    // Common Corner, analogue,no energy cut
    // createClusters
    auto mergedCells3 = Acts::createClusters<Acts::DigitizationCell>(
        testCells, nBins0, nBins1, true, true, 0.);
    // check number of clusters
    BOOST_CHECK_EQUAL(mergedCells3.size(), nClusters);

    float               data3 = 0;
    std::vector<size_t> clusterSizes3;
    for (size_t i = 0; i < mergedCells3.size(); i++) {
      auto cells = mergedCells3.at(i);
      for (auto& cell : cells) {
        data3 += cell.data;
      }
      // check the cluster sizes
      clusterSizes3.push_back(cells.size());
    }
    // check the cluster sizes
    std::sort(clusterSizes3.begin(), clusterSizes3.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizes.begin(),
                                  clusterSizes.end(),
                                  clusterSizes3.begin(),
                                  clusterSizes3.end());
    // check cells
    BOOST_CHECK_CLOSE(data3, nCells, 10e-5);

    // Common Corner, analogue, energy cut
    // createClusters
    auto mergedCells4 = Acts::createClusters<Acts::DigitizationCell>(
        testCells, nBins0, nBins1, true, true, 1.5);
    // check number of clusters
    BOOST_CHECK_EQUAL(mergedCells4.size(), nClustersCut);

    float               data4 = 0;
    std::vector<size_t> clusterSizes4;
    for (size_t i = 0; i < mergedCells4.size(); i++) {
      auto cells = mergedCells4.at(i);
      for (auto& cell : cells) {
        data4 += cell.data;
      }
      // check the cluster sizes
      clusterSizes4.push_back(cells.size());
    }
    // check the cluster sizes
    std::sort(clusterSizes4.begin(), clusterSizes4.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizesCut.begin(),
                                  clusterSizesCut.end(),
                                  clusterSizes4.begin(),
                                  clusterSizes4.end());
    // check cells
    BOOST_CHECK_CLOSE(data4, 2. * (nCells - nCellsWithoutDuplicates), 10e-5);
  }

  /// This test tests the clusterization of cells which belong to the same
  /// cluster for 8-cell/4-cell merging, digital/analogue readout and with a
  /// possible energy cut applied and a bigger grid than in create_Clusters1
  BOOST_AUTO_TEST_CASE(create_Clusters2)
  {
    std::vector<Acts::DigitizationCell> testCells1;

    size_t nCells           = 99;
    size_t delta            = 3;
    size_t nClustersNoTouch = (nCells / delta) * (nCells / delta);
    size_t nCellsInClusters = nClustersNoTouch * 3;

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
    // copy
    std::vector<Acts::DigitizationCell> testCells2 = testCells1;
    std::vector<Acts::DigitizationCell> testCells3 = testCells1;
    std::vector<Acts::DigitizationCell> testCells4 = testCells1;

    // in a 99x99 grid we get 33x33=1089 clusters
    // Now we should have the same number of cluster for common corner and
    // common edge case
    // common edge
    auto mergedCells1 = Acts::createClusters<Acts::DigitizationCell>(
        testCells1, nCells, nCells, true, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells1.size(), nClustersNoTouch);
    // common corner
    auto mergedCells2 = Acts::createClusters<Acts::DigitizationCell>(
        testCells1, nCells, nCells, false, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells2.size(), nClustersNoTouch);

    // now test merging - there is no merging at the moment
    float data1 = 0;
    for (auto& cells : mergedCells1) {
      for (auto& i : cells) {
        data1 += i.data;
        // check cluster sizes
        BOOST_CHECK_EQUAL(cells.size(), delta);
      }
    }
    BOOST_CHECK_CLOSE(data1, nCellsInClusters, 10e-5);

    // now test merging - there is no merging at the moment
    float data2 = 0;
    for (auto& cells : mergedCells2) {
      for (auto& i : cells) {
        data2 += i.data;
        // check cluster sizes
        BOOST_CHECK_EQUAL(cells.size(), delta);
      }
    }
    BOOST_CHECK_CLOSE(data2, nCellsInClusters, 10e-5);

    // now we add some cells which lead to merging only for common corner and
    // create new clusters for edge case)
    size_t delta2           = 9;
    size_t nCornerCells     = nCells / delta2 * nCells / delta2;
    size_t nClusters_merged = nClustersNoTouch - nCornerCells * 2;

    for (size_t i = 2; i < nCells; i += delta2) {
      for (size_t j = 2; j < nCells; j += delta2) {
        testCells2.push_back(Acts::DigitizationCell(i, j, 1));
      }
    }
    // common corner
    auto mergedCells3 = Acts::createClusters<Acts::DigitizationCell>(
        testCells2, nCells, nCells, true, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells3.size(), nClusters_merged);

    // common edge
    auto mergedCells4 = Acts::createClusters<Acts::DigitizationCell>(
        testCells2, nCells, nCells, false, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells4.size(), nClustersNoTouch + nCornerCells);

    // now we add some cells which lead to merging also for edge case
    for (size_t i = 2; i < nCells; i += delta2) {
      for (size_t j = 2; j < nCells; j += delta2) {
        testCells3.push_back(Acts::DigitizationCell(i, j - 1, 1));
        testCells3.push_back(Acts::DigitizationCell(i + 1, j, 1));
      }
    }

    // common corner
    auto mergedCells5 = Acts::createClusters<Acts::DigitizationCell>(
        testCells3, nCells, nCells, true, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells5.size(), nClusters_merged);

    // common edge
    auto mergedCells6 = Acts::createClusters<Acts::DigitizationCell>(
        testCells3, nCells, nCells, false, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells6.size(), nClusters_merged);

    // now adding the same cells again on two positions of clusters
    for (size_t i = 0; i < nCells; i += delta) {
      for (size_t j = 0; j < nCells; j += delta) {
        testCells4.push_back(Acts::DigitizationCell(i, j, 1));
        testCells4.push_back(Acts::DigitizationCell(i, j + 1, 1));
      }
    }

    size_t nCellsInClustersDuplicated
        = nClustersNoTouch * 3 + nClustersNoTouch * 2;

    // digital readout
    auto mergedCells7 = Acts::createClusters<Acts::DigitizationCell>(
        testCells4, nCells, nCells, true, false, 0.);
    BOOST_CHECK_EQUAL(mergedCells7.size(), nClustersNoTouch);
    float data7 = 0;
    for (auto& cells : mergedCells7) {
      for (auto& i : cells) {
        data7 += i.data;
      }
      BOOST_CHECK_EQUAL(cells.size(), 3);
    }
    BOOST_CHECK_CLOSE(data7, nCellsInClusters, 10e-5);

    // analougue readout
    auto mergedCells8 = Acts::createClusters<Acts::DigitizationCell>(
        testCells4, nCells, nCells, true, true, 0.);
    BOOST_CHECK_EQUAL(mergedCells8.size(), nClustersNoTouch);
    float data8 = 0;
    for (auto& cells : mergedCells8) {
      for (auto& i : cells) {
        data8 += i.data;
      }
      BOOST_CHECK_EQUAL(cells.size(), 3);
    }
    BOOST_CHECK_CLOSE(data8, nCellsInClustersDuplicated, 10e-5);

    // analougue readout & energy cut
    auto mergedCells9 = Acts::createClusters<Acts::DigitizationCell>(
        testCells4, nCells, nCells, true, true, 1.5);
    BOOST_CHECK_EQUAL(mergedCells9.size(), nClustersNoTouch);
    float data9 = 0;
    for (auto& cells : mergedCells9) {
      for (auto& i : cells) {
        data9 += i.data;
      }
      BOOST_CHECK_EQUAL(cells.size(), 2);
    }
    BOOST_CHECK_CLOSE(data9, (nClustersNoTouch * 2) * 2, 10e-5);
  }
}
}
