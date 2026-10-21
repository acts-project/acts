// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Clusterization/TimedClusterization.hpp"

#include <algorithm>

using namespace Acts;

// Define objects
using Identifier = std::size_t;
struct Cell {
  Cell(Identifier identifier, int c, int r, double t)
      : id(identifier), column(c), row(r), time(t) {}

  Identifier id{};
  int column{0};
  int row{0};
  int label{-1};
  double time{0.};
};

struct Cluster {
  std::vector<Identifier> ids{};
};

using CellCollection = std::vector<Cell>;
using ClusterCollection = std::vector<Cluster>;

// Define functions
static inline int getCellRow(const Cell& cell) {
  return cell.row;
}

static inline int getCellColumn(const Cell& cell) {
  return cell.column;
}

static inline double getCellTime(const Cell& cell) {
  return cell.time;
}

static void clusterAddCell(Cluster& cl, const Cell& cell) {
  cl.ids.push_back(cell.id);
}

BOOST_AUTO_TEST_SUITE(ClusterizationSuite)

BOOST_AUTO_TEST_CASE(TimedGrid_1D_withtime) {
  // 1x10 matrix
  /*
    X X X Y O X Y Y X X
  */
  // 6 + 3 cells -> 3 + 2 clusters in total

  std::vector<Cell> cells;
  // X
  cells.emplace_back(0ul, 0, -1, 0);
  cells.emplace_back(1ul, 1, -1, 0);
  cells.emplace_back(2ul, 2, -1, 0);
  cells.emplace_back(3ul, 5, -1, 0);
  cells.emplace_back(4ul, 8, -1, 0);
  cells.emplace_back(5ul, 9, -1, 0);
  // Y
  cells.emplace_back(6ul, 3, 0, 1);
  cells.emplace_back(7ul, 6, 1, 1);
  cells.emplace_back(8ul, 7, 1, 1);

  std::vector<std::vector<Identifier>> expectedResults;
  expectedResults.push_back({0ul, 1ul, 2ul});
  expectedResults.push_back({6ul});
  expectedResults.push_back({3ul});
  expectedResults.push_back({7ul, 8ul});
  expectedResults.push_back({4ul, 5ul});

  Acts::Ccl::ClusteringData data;
  ClusterCollection clusters;
  Acts::Ccl::createClusters<CellCollection, ClusterCollection, 1>(
      data, cells, clusters, Acts::Ccl::TimedConnect<Cell, 1>(0.5));

  BOOST_CHECK_EQUAL(5ul, clusters.size());

  for (std::size_t i(0); i < clusters.size(); ++i) {
    std::vector<Identifier>& timedIds = clusters[i].ids;
    const std::vector<Identifier>& expected = expectedResults[i];
    std::ranges::sort(timedIds);
    BOOST_CHECK_EQUAL(timedIds.size(), expected.size());

    for (std::size_t j(0); j < timedIds.size(); ++j) {
      BOOST_CHECK_EQUAL(timedIds[j], expected[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TimedGrid_1D_duplicate_cells) {
  // Cells with same position but time difference larger than tolerance should
  // not be considered duplicates
  CellCollection cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 0, 0, 1.0)};
  ClusterCollection clusters;
  Ccl::ClusteringData data;

  Ccl::createClusters<CellCollection, ClusterCollection, 1>(
      data, cells, clusters, Ccl::TimedConnect<Cell, 1>(0.5));
  BOOST_CHECK_EQUAL(2ul, clusters.size());

  // Cells with same position and time difference within tolerance should be
  // considered duplicates
  cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 0, 0, 0.4)};
  clusters.clear();
  data.clear();

  BOOST_CHECK_THROW(
      (Ccl::createClusters<CellCollection, ClusterCollection, 1>(
          data, cells, clusters, Ccl::TimedConnect<Cell, 1>(0.5))),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TimedGrid_1D_space_and_time) {
  // Cells with the same position but large time differences are not duplicates.
  // However, they can still be clustered through a neighboring cell if it is
  // within the time tolerance for both.
  CellCollection cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 1, 0, 0.4),
                          Cell(2ul, 0, 0, 0.8)};
  ClusterCollection clusters;
  Ccl::ClusteringData data;

  Ccl::createClusters<CellCollection, ClusterCollection, 1>(
      data, cells, clusters, Ccl::TimedConnect<Cell, 1>(0.5));
  BOOST_CHECK_EQUAL(1ul, clusters.size());
  BOOST_CHECK_EQUAL(clusters[0].ids.size(), 3ul);
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_notime) {
  // 4x4 matrix
  /*
    X O O X
    O O O X
    X X O O
    X O O X
  */
  // 7 cells -> 4 clusters in total

  std::vector<Cell> cells;
  cells.emplace_back(0ul, 0, 0, 0);
  cells.emplace_back(1ul, 3, 0, 0);
  cells.emplace_back(2ul, 3, 1, 0);
  cells.emplace_back(3ul, 0, 2, 0);
  cells.emplace_back(4ul, 1, 2, 0);
  cells.emplace_back(5ul, 0, 3, 0);
  cells.emplace_back(6ul, 3, 3, 0);

  std::vector<std::vector<Identifier>> expectedResults;
  expectedResults.push_back({0ul});
  expectedResults.push_back({3ul, 4ul, 5ul});
  expectedResults.push_back({1ul, 2ul});
  expectedResults.push_back({6ul});

  Acts::Ccl::ClusteringData data;
  ClusterCollection clusters;
  Acts::Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters,
      Acts::Ccl::TimedConnect<Cell, 2>(std::numeric_limits<double>::max()));

  BOOST_CHECK_EQUAL(4ul, clusters.size());

  // Compare against default connect (only space)
  data.clear();
  ClusterCollection defaultClusters;
  Acts::Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, defaultClusters, Acts::Ccl::DefaultConnect<Cell, 2>());

  BOOST_CHECK_EQUAL(4ul, defaultClusters.size());
  BOOST_CHECK_EQUAL(defaultClusters.size(), expectedResults.size());

  std::vector<std::size_t> sizes{1, 3, 2, 1};
  for (std::size_t i(0); i < clusters.size(); ++i) {
    std::vector<Identifier>& timedIds = clusters[i].ids;
    std::vector<Identifier>& defaultIds = defaultClusters[i].ids;
    const std::vector<Identifier>& expected = expectedResults[i];
    BOOST_CHECK_EQUAL(timedIds.size(), defaultIds.size());
    BOOST_CHECK_EQUAL(timedIds.size(), sizes[i]);
    BOOST_CHECK_EQUAL(timedIds.size(), expected.size());

    std::ranges::sort(timedIds);
    std::ranges::sort(defaultIds);
    for (std::size_t j(0); j < timedIds.size(); ++j) {
      BOOST_CHECK_EQUAL(timedIds[j], defaultIds[j]);
      BOOST_CHECK_EQUAL(timedIds[j], expected[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_withtime) {
  // 4x4 matrix
  /*
    X Y O X
    O Y Y X
    X X Z Z
    X O O X
  */
  // 7 + 3 + 2 cells -> 4 + 1 + 1 clusters in total

  std::vector<Cell> cells;
  // X
  cells.emplace_back(0ul, 0, 0, 0);
  cells.emplace_back(1ul, 3, 0, 0);
  cells.emplace_back(2ul, 3, 1, 0);
  cells.emplace_back(3ul, 0, 2, 0);
  cells.emplace_back(4ul, 1, 2, 0);
  cells.emplace_back(5ul, 0, 3, 0);
  cells.emplace_back(6ul, 3, 3, 0);
  // Y
  cells.emplace_back(7ul, 1, 0, 1);
  cells.emplace_back(8ul, 1, 1, 1);
  cells.emplace_back(9ul, 2, 1, 1);
  // Z
  cells.emplace_back(10ul, 2, 2, 2);
  cells.emplace_back(11ul, 3, 2, 2);

  std::vector<std::vector<Identifier>> expectedResults;
  expectedResults.push_back({0ul});
  expectedResults.push_back({3ul, 4ul, 5ul});
  expectedResults.push_back({7ul, 8ul, 9ul});
  expectedResults.push_back({10ul, 11ul});
  expectedResults.push_back({1ul, 2ul});
  expectedResults.push_back({6ul});

  Acts::Ccl::ClusteringData data;
  ClusterCollection clusters;
  Acts::Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters, Acts::Ccl::TimedConnect<Cell, 2>(0.5));

  BOOST_CHECK_EQUAL(6ul, clusters.size());

  std::vector<std::size_t> sizes{1, 3, 3, 2, 2, 1};
  for (std::size_t i(0); i < clusters.size(); ++i) {
    std::vector<Identifier>& timedIds = clusters[i].ids;
    BOOST_CHECK_EQUAL(timedIds.size(), sizes[i]);
    std::ranges::sort(timedIds);

    const std::vector<Identifier>& expected = expectedResults[i];
    BOOST_CHECK_EQUAL(timedIds.size(), expected.size());

    for (std::size_t j(0); j < timedIds.size(); ++j) {
      BOOST_CHECK_EQUAL(timedIds[j], expected[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_noTollerance) {
  // 4x4 matrix
  /*
    X O O X
    O O O X
    X X O O
    X O O X
   */
  // 7 cells -> 7 clusters in total
  // since time requirement will never be satisfied

  std::vector<Cell> cells;
  cells.emplace_back(0ul, 0, 0, 0);
  cells.emplace_back(1ul, 3, 0, 0);
  cells.emplace_back(2ul, 3, 1, 0);
  cells.emplace_back(3ul, 0, 2, 0);
  cells.emplace_back(4ul, 1, 2, 0);
  cells.emplace_back(5ul, 0, 3, 0);
  cells.emplace_back(6ul, 3, 3, 0);

  std::vector<std::vector<Identifier>> expectedResults;
  expectedResults.push_back({0ul});
  expectedResults.push_back({3ul});
  expectedResults.push_back({5ul});
  expectedResults.push_back({4ul});
  expectedResults.push_back({1ul});
  expectedResults.push_back({2ul});
  expectedResults.push_back({6ul});

  Acts::Ccl::ClusteringData data;
  ClusterCollection clusters;
  Acts::Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters, Acts::Ccl::TimedConnect<Cell, 2>(0.));

  BOOST_CHECK_EQUAL(7ul, clusters.size());

  for (std::size_t i(0); i < clusters.size(); ++i) {
    std::vector<Identifier>& timedIds = clusters[i].ids;
    const std::vector<Identifier>& expected = expectedResults[i];

    BOOST_CHECK_EQUAL(timedIds.size(), 1);
    BOOST_CHECK_EQUAL(timedIds.size(), expected.size());
    BOOST_CHECK_EQUAL(timedIds[0], expected[0]);
  }
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_duplicate_cells) {
  // Cells with same position but time difference larger than tolerance should
  // not be considered duplicates
  CellCollection cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 0, 0, 1.0)};
  ClusterCollection clusters;
  Ccl::ClusteringData data;

  Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters, Ccl::TimedConnect<Cell, 2>(0.5));
  BOOST_CHECK_EQUAL(2ul, clusters.size());

  // Cells with same position and time difference within tolerance should be
  // considered duplicates
  cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 0, 0, 0.4)};
  clusters.clear();
  data.clear();

  BOOST_CHECK_THROW(
      (Ccl::createClusters<CellCollection, ClusterCollection, 2>(
          data, cells, clusters, Ccl::TimedConnect<Cell, 2>(0.5))),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_space_and_time_conn4) {
  // Cells with the same position but large time differences are not duplicates.
  // However, they can still be clustered through a neighboring cell if it is
  // within the time tolerance for both:
  // 3x3 matrix
  /*
    O  X/Z  O
   X/Z  Y  X/Z
    O  X/Z  O
  */
  CellCollection cells = {
      Cell(0ul, 1, 0, 0.0), Cell(1ul, 0, 1, 0.0), Cell(2ul, 2, 1, 0.0),
      Cell(3ul, 1, 2, 0.0), Cell(4ul, 1, 0, 0.8), Cell(5ul, 0, 1, 0.8),
      Cell(6ul, 2, 1, 0.8), Cell(7ul, 1, 2, 0.8), Cell(8ul, 1, 1, 0.4)};
  ClusterCollection clusters;
  Ccl::ClusteringData data;

  Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters, Ccl::TimedConnect<Cell, 2>(0.5, false));
  BOOST_CHECK_EQUAL(1ul, clusters.size());
  BOOST_CHECK_EQUAL(clusters[0].ids.size(), 9ul);
}

BOOST_AUTO_TEST_CASE(TimedGrid_2D_space_and_time_conn8) {
  // Cells with the same position but large time differences are not duplicates.
  // However, they can still be clustered through a neighboring cell if it is
  // within the time tolerance for both:
  // 3x3 matrix
  /*
    Z   X   O
    X   Y   O
   X/Z  O   O
  */
  CellCollection cells = {Cell(0ul, 0, 0, 0.0), Cell(1ul, 0, 2, 0.0),
                          Cell(2ul, 1, 0, 0.8), Cell(3ul, 0, 1, 0.8),
                          Cell(4ul, 0, 2, 0.8), Cell(5ul, 1, 1, 0.4)};
  ClusterCollection clusters;
  Ccl::ClusteringData data;

  Ccl::createClusters<CellCollection, ClusterCollection, 2>(
      data, cells, clusters, Ccl::TimedConnect<Cell, 2>(0.5, true));
  BOOST_CHECK_EQUAL(1ul, clusters.size());
  BOOST_CHECK_EQUAL(clusters[0].ids.size(), 6ul);
}

BOOST_AUTO_TEST_SUITE_END()
