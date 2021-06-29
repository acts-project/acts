// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Clusterization/Clusterization.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <utility>
#include <vector>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

struct Cell {
  Cell(size_t c0, size_t c1, float act)
      : channel0(c0), channel1(c1), data(act) {}

  void addCell(const Cell& c, bool analogue) {
    if (analogue)
      data += c.data;
  }
  double activation() const { return data; }

  size_t channel0;
  size_t channel1;
  double data;
};

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
BOOST_AUTO_TEST_CASE(create_Clusters1) {
  size_t nBins0 = 10;
  std::vector<size_t> clusterSizes = {2, 2, 2, 3, 6, 6, 7};
  std::vector<size_t> clusterSizesCut = {1, 1, 1, 1, 1, 1, 1, 2, 2, 3};
  std::vector<size_t> clusterSizesEdge = {1, 1, 1, 1, 1, 1, 1,
                                          2, 2, 2, 3, 3, 3, 6};
  size_t nClusters = clusterSizes.size();
  size_t nClustersCut = clusterSizesCut.size();
  size_t nClustersEdge = clusterSizesEdge.size();

  auto globalIndex = [&nBins0](const size_t& a, const size_t& b) {
    return (a + nBins0 * b);
  };

  std::unordered_map<size_t, std::pair<Cell, bool>> testCells;
  // add cells covering all cases
  testCells.insert({globalIndex(0, 0), {Cell(0, 0, 1), false}});
  testCells.insert({globalIndex(0, 3), {Cell(0, 3, 1), false}});
  testCells.insert({globalIndex(0, 4), {Cell(0, 4, 1), false}});
  testCells.insert({globalIndex(0, 8), {Cell(0, 8, 1), false}});
  testCells.insert({globalIndex(1, 1), {Cell(1, 1, 1), false}});
  testCells.insert({globalIndex(1, 8), {Cell(1, 8, 1), false}});
  testCells.insert({globalIndex(2, 5), {Cell(2, 5, 1), false}});
  testCells.insert({globalIndex(2, 8), {Cell(2, 8, 1), false}});
  testCells.insert({globalIndex(3, 1), {Cell(3, 1, 1), false}});
  testCells.insert({globalIndex(3, 4), {Cell(3, 4, 1), false}});
  testCells.insert({globalIndex(3, 6), {Cell(3, 6, 1), false}});
  testCells.insert({globalIndex(4, 0), {Cell(4, 0, 1), false}});
  testCells.insert({globalIndex(4, 6), {Cell(4, 6, 1), false}});
  testCells.insert({globalIndex(5, 3), {Cell(5, 3, 1), false}});
  testCells.insert({globalIndex(5, 5), {Cell(5, 5, 1), false}});
  testCells.insert({globalIndex(5, 4), {Cell(5, 4, 1), false}});
  testCells.insert({globalIndex(5, 5), {Cell(5, 5, 1), false}});
  testCells.insert({globalIndex(7, 1), {Cell(7, 1, 1), false}});
  testCells.insert({globalIndex(7, 2), {Cell(7, 2, 1), false}});
  testCells.insert({globalIndex(7, 5), {Cell(7, 5, 1), false}});
  testCells.insert({globalIndex(7, 6), {Cell(7, 6, 1), false}});
  testCells.insert({globalIndex(8, 2), {Cell(8, 2, 1), false}});
  testCells.insert({globalIndex(8, 7), {Cell(8, 7, 1), false}});
  testCells.insert({globalIndex(9, 0), {Cell(9, 0, 1), false}});
  testCells.insert({globalIndex(9, 1), {Cell(9, 1, 1), false}});
  testCells.insert({globalIndex(9, 2), {Cell(9, 2, 1), false}});
  testCells.insert({globalIndex(9, 4), {Cell(9, 4, 1), false}});
  testCells.insert({globalIndex(9, 5), {Cell(9, 5, 1), false}});
  testCells.insert({globalIndex(9, 6), {Cell(9, 6, 1), false}});

  size_t nCellsWithoutDuplicates = testCells.size();

  std::unordered_map<size_t, std::pair<Cell, bool>> testCells1;
  // add cells covering all cases
  testCells1.insert({globalIndex(0, 0), {Cell(0, 0, 1), false}});
  testCells1.insert({globalIndex(0, 3), {Cell(0, 3, 2), false}});
  testCells1.insert({globalIndex(0, 4), {Cell(0, 4, 1), false}});
  testCells1.insert({globalIndex(0, 8), {Cell(0, 8, 1), false}});
  testCells1.insert({globalIndex(1, 1), {Cell(1, 1, 2), false}});
  testCells1.insert({globalIndex(1, 8), {Cell(1, 8, 2), false}});
  testCells1.insert({globalIndex(2, 5), {Cell(2, 5, 2), false}});
  testCells1.insert({globalIndex(2, 8), {Cell(2, 8, 2), false}});
  testCells1.insert({globalIndex(3, 1), {Cell(3, 1, 1), false}});
  testCells1.insert({globalIndex(3, 4), {Cell(3, 4, 1), false}});
  testCells1.insert({globalIndex(3, 6), {Cell(3, 6, 1), false}});
  testCells1.insert({globalIndex(4, 0), {Cell(4, 0, 2), false}});
  testCells1.insert({globalIndex(4, 6), {Cell(4, 6, 1), false}});
  testCells1.insert({globalIndex(5, 3), {Cell(5, 3, 1), false}});
  testCells1.insert({globalIndex(5, 5), {Cell(5, 5, 1), false}});
  testCells1.insert({globalIndex(5, 4), {Cell(5, 4, 2), false}});
  testCells1.insert({globalIndex(5, 5), {Cell(5, 5, 1), false}});
  testCells1.insert({globalIndex(7, 1), {Cell(7, 1, 2), false}});
  testCells1.insert({globalIndex(7, 2), {Cell(7, 2, 1), false}});
  testCells1.insert({globalIndex(7, 5), {Cell(7, 5, 1), false}});
  testCells1.insert({globalIndex(7, 6), {Cell(7, 6, 2), false}});
  testCells1.insert({globalIndex(8, 2), {Cell(8, 2, 1), false}});
  testCells1.insert({globalIndex(8, 7), {Cell(8, 7, 2), false}});
  testCells1.insert({globalIndex(9, 0), {Cell(9, 0, 2), false}});
  testCells1.insert({globalIndex(9, 1), {Cell(9, 1, 2), false}});
  testCells1.insert({globalIndex(9, 2), {Cell(9, 2, 2), false}});
  testCells1.insert({globalIndex(9, 4), {Cell(9, 4, 2), false}});
  testCells1.insert({globalIndex(9, 5), {Cell(9, 5, 1), false}});
  testCells1.insert({globalIndex(9, 6), {Cell(9, 6, 1), false}});

  // add duplicates

  size_t nCells = testCells1.size();

  // make copies for all the tests (the boolean in the map gets changed)
  auto testCells2 = testCells;
  auto testCells3 = testCells;
  auto testCells4 = testCells1;

  // Common Corner, digital,no energy cut
  // createClusters
  auto mergedCells1 = Acts::createClusters<Cell>(testCells, nBins0, true, 0.);
  // check number of clusters
  BOOST_CHECK_EQUAL(mergedCells1.size(), nClusters);

  float data1 = 0;
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
  BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizes.begin(), clusterSizes.end(),
                                clusterSizes1.begin(), clusterSizes1.end());
  // check cells
  CHECK_CLOSE_REL(data1, nCellsWithoutDuplicates, 1e-5);

  // Common Edge, digital,no energy cut
  // createClusters
  auto mergedCells2 = Acts::createClusters<Cell>(testCells2, nBins0, false, 0.);
  // check number of clusters
  BOOST_CHECK_EQUAL(mergedCells2.size(), nClustersEdge);

  float data2 = 0;
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
                                clusterSizesEdge.end(), clusterSizes2.begin(),
                                clusterSizes2.end());
  // check cells
  CHECK_CLOSE_REL(data2, nCellsWithoutDuplicates, 1e-5);

  // Common Corner, analogue,no energy cut
  // createClusters
  auto mergedCells3 = Acts::createClusters<Cell>(testCells3, nBins0, true, 0.);
  // check number of clusters
  BOOST_CHECK_EQUAL(mergedCells3.size(), nClusters);

  float data3 = 0;
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
  BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizes.begin(), clusterSizes.end(),
                                clusterSizes3.begin(), clusterSizes3.end());
  // check cells
  CHECK_CLOSE_REL(data3, nCells, 1e-5);

  // Common Corner, analogue, energy cut
  // createClusters
  auto mergedCells4 = Acts::createClusters<Cell>(testCells4, nBins0, true, 1.5);
  // check number of clusters

  BOOST_CHECK_EQUAL(mergedCells4.size(), nClustersCut);

  float data4 = 0;
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
  BOOST_CHECK_EQUAL_COLLECTIONS(clusterSizesCut.begin(), clusterSizesCut.end(),
                                clusterSizes4.begin(), clusterSizes4.end());
  // check cells
  CHECK_CLOSE_REL(data4, nCells, 1e-5);
}

/// This test tests the clusterization of cells which belong to the same
/// cluster for 8-cell/4-cell merging, digital/analogue readout and with a
/// possible energy cut applied and a bigger grid than in create_Clusters1
BOOST_AUTO_TEST_CASE(create_Clusters2) {
  std::unordered_map<size_t, std::pair<Cell, bool>> testCells1;

  size_t nCells = 99;
  size_t delta = 3;
  size_t nClustersNoTouch = (nCells / delta) * (nCells / delta);
  size_t nCellsInClusters = nClustersNoTouch * 3;

  auto globalIndex = [&nCells](const size_t& a, const size_t& b) {
    return (a + nCells * b);
  };

  // Create clusters which are separated by one cell always
  // Clusters have the following shape:
  // -----
  // --##-
  // ---#-
  // -----
  for (size_t i = 0; i < nCells; i += delta) {
    for (size_t j = 0; j < nCells; j += delta) {
      auto cellA = Cell(i, j, 1);
      auto cellB = Cell(i, j + 1, 1);
      auto cellC = Cell(i + 1, j + 1, 1);

      auto insertCellA = testCells1.insert({globalIndex(i, j), {cellA, false}});
      if (!insertCellA.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellA.first->second.first.addCell(cellA, false);
      }
      auto insertCellB =
          testCells1.insert({globalIndex(i, j + 1), {cellB, false}});
      if (!insertCellB.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellB.first->second.first.addCell(cellB, false);
      }
      auto insertCellC =
          testCells1.insert({globalIndex(i + 1, j + 1), {cellC, false}});
      if (!insertCellC.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellC.first->second.first.addCell(cellC, false);
      }
    }
  }
  // copy
  auto testCells2 = testCells1;
  auto testCells3 = testCells1;
  auto testCells5 = testCells1;
  auto testCells7 = testCells1;
  auto testCells8 = testCells1;

  // in a 99x99 grid we get 33x33=1089 clusters
  // Now we should have the same number of cluster for common corner and
  // common edge case
  // common edge
  auto mergedCells1 = Acts::createClusters<Cell>(testCells1, nCells, true, 0.);
  BOOST_CHECK_EQUAL(mergedCells1.size(), nClustersNoTouch);
  // common corner
  auto mergedCells2 = Acts::createClusters<Cell>(testCells2, nCells, false, 0.);
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
  CHECK_CLOSE_REL(data1, nCellsInClusters, 1e-5);

  // now test merging - there is no merging at the moment
  float data2 = 0;
  for (auto& cells : mergedCells2) {
    for (auto& i : cells) {
      data2 += i.data;
      // check cluster sizes
      BOOST_CHECK_EQUAL(cells.size(), delta);
    }
  }
  CHECK_CLOSE_REL(data2, nCellsInClusters, 1e-5);

  // now we add some cells which lead to merging only for common corner and
  // create new clusters for edge case)
  size_t delta2 = 9;
  size_t nCornerCells = nCells / delta2 * nCells / delta2;
  size_t nClusters_merged = nClustersNoTouch - nCornerCells * 2;

  for (size_t i = 2; i < nCells; i += delta2) {
    for (size_t j = 2; j < nCells; j += delta2) {
      auto cell = Cell(i, j, 1);
      auto insertCell = testCells3.insert({globalIndex(i, j), {cell, false}});
      if (!insertCell.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCell.first->second.first.addCell(cell, false);
      }
    }
  }
  auto testCells4 = testCells3;

  // common corner
  auto mergedCells3 = Acts::createClusters<Cell>(testCells3, nCells, true, 0.);
  BOOST_CHECK_EQUAL(mergedCells3.size(), nClusters_merged);

  // common edge
  auto mergedCells4 = Acts::createClusters<Cell>(testCells4, nCells, false, 0.);
  BOOST_CHECK_EQUAL(mergedCells4.size(), nClustersNoTouch + nCornerCells);

  // now we add some cells which lead to merging also for edge case
  for (size_t i = 2; i < nCells; i += delta2) {
    for (size_t j = 2; j < nCells; j += delta2) {
      auto cellA = Cell(i, j - 1, 1);
      auto insertCellA =
          testCells5.insert({globalIndex(i, j - 1), {cellA, false}});
      if (!insertCellA.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellA.first->second.first.addCell(cellA, false);
      }
      auto cellB = Cell(i + 1, j, 1);
      auto insertCellB =
          testCells5.insert({globalIndex(i + 1, j), {cellB, false}});
      if (!insertCellB.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellB.first->second.first.addCell(cellB, false);
      }
    }
  }
  auto testCells6 = testCells5;

  // common corner
  auto mergedCells5 = Acts::createClusters<Cell>(testCells5, nCells, true, 0.);
  BOOST_CHECK_EQUAL(mergedCells5.size(), nClusters_merged);

  // common edge
  auto mergedCells6 = Acts::createClusters<Cell>(testCells6, nCells, false, 0.);
  BOOST_CHECK_EQUAL(mergedCells6.size(), nClusters_merged);

  // now adding the same cells again on two positions of clusters - digital
  // readout
  for (size_t i = 0; i < nCells; i += delta) {
    for (size_t j = 0; j < nCells; j += delta) {
      auto cellA = Cell(i, j, 1);
      auto cellB = Cell(i, j + 1, 1);

      auto insertCellA = testCells7.insert({globalIndex(i, j), {cellA, false}});
      if (!insertCellA.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellA.first->second.first.addCell(cellA, false);
      }
      auto insertCellB =
          testCells7.insert({globalIndex(i, j + 1), {cellB, false}});
      if (!insertCellB.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellB.first->second.first.addCell(cellB, false);
      }
      // analogue readout
      insertCellA = testCells8.insert({globalIndex(i, j), {cellA, false}});
      if (!insertCellA.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellA.first->second.first.addCell(cellA, true);
      }
      insertCellB = testCells8.insert({globalIndex(i, j + 1), {cellB, false}});
      if (!insertCellB.second) {
        // check if there is already a cell at same position and merge in that
        // case
        insertCellB.first->second.first.addCell(cellB, true);
      }
    }
  }
  auto testCells9 = testCells8;

  size_t nCellsInClustersDuplicated =
      nClustersNoTouch * 3 + nClustersNoTouch * 2;

  // digital readout
  auto mergedCells7 = Acts::createClusters<Cell>(testCells7, nCells, true, 0.);
  BOOST_CHECK_EQUAL(mergedCells7.size(), nClustersNoTouch);
  float data7 = 0;
  for (auto& cells : mergedCells7) {
    for (auto& i : cells) {
      data7 += i.data;
    }
    BOOST_CHECK_EQUAL(cells.size(), 3u);
  }
  CHECK_CLOSE_REL(data7, nCellsInClusters, 1e-5);

  // analougue readout
  auto mergedCells8 = Acts::createClusters<Cell>(testCells8, nCells, true, 0.);
  BOOST_CHECK_EQUAL(mergedCells8.size(), nClustersNoTouch);
  float data8 = 0;
  for (auto& cells : mergedCells8) {
    for (auto& i : cells) {
      data8 += i.data;
    }
    BOOST_CHECK_EQUAL(cells.size(), 3u);
  }
  CHECK_CLOSE_REL(data8, nCellsInClustersDuplicated, 1e-5);

  // analougue readout & energy cut
  auto mergedCells9 = Acts::createClusters<Cell>(testCells9, nCells, true, 1.5);
  BOOST_CHECK_EQUAL(mergedCells9.size(), nClustersNoTouch);
  float data9 = 0;
  for (auto& cells : mergedCells9) {
    for (auto& i : cells) {
      data9 += i.data;
    }
    BOOST_CHECK_EQUAL(cells.size(), 2u);
  }
  CHECK_CLOSE_REL(data9, (nClustersNoTouch * 2) * 2, 1e-5);
}
}  // namespace Test
}  // namespace Acts
