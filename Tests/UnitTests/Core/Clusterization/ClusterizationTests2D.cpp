// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include <boost/test/unit_test.hpp>

#include "Acts/Clusterization/Clusterization.hpp"
#include "Acts/Clusterization/InPlaceClusterization.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

#include "ClusterizationTestData.hpp"
#include "ClusterizationTestData8Fold.hpp"

using namespace Acts;

using Rectangle = std::array<int, 4>;

std::vector<Rectangle> concat(std::vector<std::vector<Rectangle>> vecs) {
  std::vector<Rectangle> flat;
  for (std::vector<Rectangle>& v : vecs) {
    if (flat.empty()) {
      flat = std::move(v);
    } else {
      flat.reserve(flat.size() + v.size());
      std::move(v.begin(), v.end(), std::back_inserter(flat));
      v.clear();
    }
  }
  return flat;
}

template <typename RNG>
std::vector<Rectangle> segment(int x0, int y0, int x1, int y1, RNG& rng) {
  // Leave enough space for one-cell buffer around clusters
  int xmin = x0 + 3;
  int xmax = x1 - 3;
  int ymin = y0 + 3;
  int ymax = y1 - 3;

  // terminal case 1
  if (xmax < xmin || ymax < ymin) {
    return {{x0, y0, x1, y1}};
  }

  std::bernoulli_distribution cointoss;
  bool splitx = cointoss(rng);
  bool splity = cointoss(rng);

  // terminal case 2
  if (!(splitx || splity)) {
    return {{x0, y0, x1, y1}};
  }

  int x_ = std::uniform_int_distribution<std::int32_t>(xmin, xmax)(rng);
  int y_ = std::uniform_int_distribution<std::int32_t>(ymin, ymax)(rng);

  if (splitx && !splity) {
    return concat({segment(x0, y0, x_, y1, rng), segment(x_, y0, x1, y1, rng)});
  } else if (!splitx && splity) {
    return concat({segment(x0, y0, x1, y_, rng), segment(x0, y_, x1, y1, rng)});
  } else if (splitx && splity) {
    return concat({segment(x0, y0, x_, y_, rng), segment(x_, y0, x1, y_, rng),
                   segment(x0, y_, x_, y1, rng), segment(x_, y_, x1, y1, rng)});
  }
  throw std::runtime_error("unreachable");
}

struct Cell2D {
  Cell2D(int rowv, int colv) : row(rowv), col(colv) {}
  int row, col;
  Ccl::Label label{Ccl::NO_LABEL};
};

int getCellRow(const Cell2D& cell) {
  return cell.row;
}

int getCellColumn(const Cell2D& cell) {
  return cell.col;
}

bool operator==(const Cell2D& left, const Cell2D& right) {
  return left.row == right.row && left.col == right.col;
}

bool cellComp(const Cell2D& left, const Cell2D& right) {
  return (left.row == right.row) ? left.col < right.col : left.row < right.row;
}

struct Cluster2D {
  std::vector<Cell2D> cells;
  std::size_t hash{0};
};

void clusterAddCell(Cluster2D& cl, const Cell2D& cell) {
  cl.cells.push_back(cell);
}

void hash(Cluster2D& cl) {
  std::ranges::sort(cl.cells, cellComp);
  cl.hash = 0;
  for (const Cell2D& c : cl.cells) {
    boost::hash_combine(cl.hash, c.col);
  }
}

bool clHashComp(const Cluster2D& left, const Cluster2D& right) {
  return left.hash < right.hash;
}

template <typename RNG>
void genclusterw(int x, int y, int x0, int y0, int x1, int y1,
                 std::vector<Cell2D>& cells, RNG& rng, double startp = 0.5,
                 double decayp = 0.9) {
  std::vector<Cell2D> add;

  auto maybe_add = [&](int x_, int y_) {
    Cell2D c(x_, y_);
    if (std::uniform_real_distribution<double>()(rng) < startp &&
        !rangeContainsValue(cells, c)) {
      cells.push_back(c);
      add.push_back(c);
    }
  };

  // NORTH
  if (y < y1) {
    maybe_add(x, y + 1);
  }
  // NORTHEAST
  if (x < x1 && y < y1) {
    maybe_add(x + 1, y + 1);
  }
  // EAST
  if (x < x1) {
    maybe_add(x + 1, y);
  }
  // SOUTHEAST
  if (x < x1 && y > y0) {
    maybe_add(x + 1, y - 1);
  }
  // SOUTH
  if (y > y0) {
    maybe_add(x, y - 1);
  }
  // SOUTHWEST
  if (x > x0 && y > y0) {
    maybe_add(x - 1, y - 1);
  }
  // WEST
  if (x > x0) {
    maybe_add(x - 1, y);
  }
  // NORTHWEST
  if (x > x0 && y < y1) {
    maybe_add(x - 1, y + 1);
  }

  for (Cell2D& c : add) {
    genclusterw(c.row, c.col, x0, y0, x1, y1, cells, rng, startp * decayp,
                decayp);
  }
}

template <typename RNG>
Cluster2D gencluster(int x0, int y0, int x1, int y1, RNG& rng,
                     double startp = 0.5, double decayp = 0.9) {
  int x0_ = x0 + 1;
  int x1_ = x1 - 1;
  int y0_ = y0 + 1;
  int y1_ = y1 - 1;

  int x = std::uniform_int_distribution<std::int32_t>(x0_, x1_)(rng);
  int y = std::uniform_int_distribution<std::int32_t>(y0_, y1_)(rng);

  std::vector<Cell2D> cells = {Cell2D(x, y)};
  genclusterw(x, y, x0_, y0_, x1_, y1_, cells, rng, startp, decayp);

  Cluster2D cl;
  cl.cells = std::move(cells);

  return cl;
}

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(ClusterizationSuite)

BOOST_AUTO_TEST_CASE(Grid_2D_rand) {
  using Cell = Cell2D;
  using CellC = std::vector<Cell>;
  using Cluster = Cluster2D;
  using ClusterC = std::vector<Cluster>;

  std::size_t sizeX = 1000;
  std::size_t sizeY = 1000;
  std::size_t startSeed = 71902647;
  std::size_t ntries = 100;

  std::cout << "Grid_2D_rand test with parameters: " << std::endl;
  std::cout << " sizeX = " << sizeX << std::endl;
  std::cout << " sizeY = " << sizeY << std::endl;
  std::cout << " startSeed = " << startSeed << std::endl;
  std::cout << " ntries = " << ntries << std::endl;

  while (ntries-- > 0) {
    std::mt19937_64 rnd(startSeed++);

    std::vector<Cluster> cls;
    std::vector<Cell> cells;
    for (Rectangle& rect : segment(0, 0, sizeX, sizeY, rnd)) {
      auto& [x0, y0, x1, y1] = rect;
      Cluster cl = gencluster(x0, y0, x1, y1, rnd);
      hash(cl);
      cells.insert(cells.end(), cl.cells.begin(), cl.cells.end());
      cls.push_back(cl);
    }

    std::shuffle(cells.begin(), cells.end(), rnd);

    Ccl::ClusteringData data;
    ClusterC newCls;
    Ccl::createClusters<CellC, ClusterC>(data, cells, newCls);

    for (Cluster& cl : newCls) {
      hash(cl);
    }

    std::ranges::sort(cls, clHashComp);
    std::ranges::sort(newCls, clHashComp);

    BOOST_CHECK_EQUAL(cls.size(), newCls.size());
    for (std::size_t i = 0; i < cls.size(); i++) {
      BOOST_CHECK_EQUAL(cls.at(i).hash, newCls.at(i).hash);
    }
  }
}

std::vector<Cluster2D> createClusterTestData(std::int16_t* test_data) {
  std::vector<Cluster2D> clusters;
  unsigned int begin_idx = 0;

  static constexpr std::int16_t end_marker = 32767;
  unsigned int n_end_marker = 0;
  for (unsigned int idx = 0; n_end_marker < 4; ++idx) {
    if (test_data[idx] == end_marker) {
      ++n_end_marker;
    } else if (n_end_marker == 2) {
      if (begin_idx + 2 == idx)
        break;
      clusters.emplace_back();
      for (unsigned int cluster_idx = begin_idx; cluster_idx < idx - 2;
           cluster_idx += 2) {
        clusters.back().cells.emplace_back(test_data[cluster_idx],
                                           test_data[cluster_idx + 1]);
      }
      std::ranges::sort(clusters.back().cells, cellComp);
      hash(clusters.back());
      begin_idx = idx;
      n_end_marker = 0;
    }
  }
  return clusters;
}

BOOST_AUTO_TEST_CASE(ClusterizeTestData) {
  // @TODO Currently fails for out-commented test_data_sets
  std::vector<std::pair<bool, std::int16_t*>> test_data_sets{
      // std::make_pair(false, cluster_test_data),
      std::make_pair(true, cluster_test_data8)
      // ,std::make_pair(true, cluster_test_data8_2)
  };
  for (auto [fold8, test_data_input] : test_data_sets) {
    using Cell = Cell2D;
    using CellC = std::vector<Cell>;
    using Cluster = Cluster2D;
    using ClusterC = std::vector<Cluster>;

    std::vector<Cluster> cls = createClusterTestData(test_data_input);
    std::vector<Cell> cells;
    for (const Cluster& a_cluster : cls) {
      cells.insert(cells.end(), a_cluster.cells.begin(), a_cluster.cells.end());
    }

    std::size_t startSeed = 71902647;
    std::mt19937_64 rnd(startSeed++);
    std::shuffle(cells.begin(), cells.end(), rnd);

    Ccl::ClusteringData data;
    ClusterC newCls;
    Ccl::createClusters<CellC, ClusterC>(
        data, cells, newCls, Acts::Ccl::DefaultConnect<Cell, 2>(fold8));

    for (Cluster& cl : newCls) {
      hash(cl);
    }

    std::ranges::sort(cls, clHashComp);
    std::ranges::sort(newCls, clHashComp);
    std::unordered_map<std::uint64_t, unsigned int> cl_map;
    cl_map.reserve(newCls.size());
    auto makeKey = [](const Cluster2D& a_cluster) {
      return (static_cast<std::uint64_t>(a_cluster.cells.front().row) << 32) |
             a_cluster.cells.front().col;
    };
    {
      unsigned int cluster_idx = 0;
      for (const Cluster2D& a_cluster : newCls) {
        cl_map.insert(std::make_pair(makeKey(a_cluster), cluster_idx));
        ++cluster_idx;
      }
    }
    BOOST_CHECK_EQUAL(cls.size(), newCls.size());
    for (std::size_t i = 0; i < cls.size(); i++) {
      std::unordered_map<std::uint64_t, unsigned int>::const_iterator iter =
          cl_map.find(makeKey(cls.at(i)));
      BOOST_CHECK_EQUAL(iter != cl_map.end(), true);
      if (iter != cl_map.end()) {
        const Cluster2D& new_cl = newCls.at(iter->second);

        BOOST_CHECK_EQUAL(cls.at(i).hash, new_cl.hash);
      }
    }
  }
}

std::vector<std::size_t> computeClusterHashes(
    std::vector<InPlaceClusterization::Cell<std::int16_t, 2, std::uint16_t>>&
        cells) {
  using T_Cell = InPlaceClusterization::Cell<std::int16_t, 2, std::uint16_t>;
  using T_CellCollection = std::vector<T_Cell>;
  std::vector<std::size_t> cluster_hashes;
  cluster_hashes.reserve(cells.size());
  for_each_cluster(
      cells, [&cluster_hashes](T_CellCollection& all_cells,
                               unsigned int idx_begin, unsigned int idx_end) {
        std::span<T_Cell> cell_range(all_cells.data() + idx_begin,
                                     all_cells.data() + idx_end);
        std::ranges::sort(cell_range, [](const T_Cell& a, const T_Cell& b) {
          return a.m_coordinates < b.m_coordinates;
        });
        std::size_t hash{};
        for (const T_Cell& a_cell : cell_range) {
          boost::hash_combine(hash, a_cell.m_coordinates[1]);
        }
        // std::cout << "Cluster " << cell_range.size() << " " << hash <<
        // std::endl;
        cluster_hashes.push_back(hash);
      });
  return cluster_hashes;
}

//@TODO also add a test for the clusterizing cells on a 1D grid.

// Test for the in-place clusterization algorithm using generated clusters
BOOST_AUTO_TEST_CASE(InPlaceClusterization_Grid_2D_rand) {
  static constexpr unsigned int SORT_AXIS = 0;
  using Cell = InPlaceClusterization::Cell<std::int16_t, 2, std::uint16_t>;
  using Cluster = Cluster2D;

  std::size_t sizeX = 1000;
  std::size_t sizeY = 1000;
  std::size_t startSeed = 71902647;
  std::size_t ntries = 100;

  std::cout << "InPlaceClusterization_Grid_2D_rand test with parameters: "
            << std::endl;
  std::cout << " sizeX = " << sizeX << std::endl;
  std::cout << " sizeY = " << sizeY << std::endl;
  std::cout << " startSeed = " << startSeed << std::endl;
  std::cout << " ntries = " << ntries << std::endl;

  while (ntries-- > 0) {
    std::mt19937_64 rnd(startSeed++);

    std::vector<Cluster> cls;
    std::vector<Cell> cells;
    unsigned int cell_idx = 0;
    for (Rectangle& rect : segment(0, 0, sizeX, sizeY, rnd)) {
      auto& [x0, y0, x1, y1] = rect;
      Cluster cl = gencluster(x0, y0, x1, y1, rnd);
      hash(cl);

      // grow cell collection with cells of generated cluster
      for (const Cell2D& a_cell : cl.cells) {
        // make sure that the coordinates in the range supported by the chosen
        // coordinate type
        assert(static_cast<std::int16_t>(getCellRow(a_cell)) ==
               getCellRow(a_cell));
        assert(static_cast<std::int16_t>(getCellColumn(a_cell)) ==
               getCellColumn(a_cell));
        assert(static_cast<std::uint16_t>(cell_idx) == cell_idx);
        cells.push_back(Cell(
            std::array<std::int16_t, 2>{
                static_cast<std::int16_t>(getCellRow(a_cell)),
                static_cast<std::int16_t>(getCellColumn(a_cell))},
            static_cast<std::uint16_t>(cell_idx)));
        ++cell_idx;
      }
      cls.push_back(cl);
    }

    // randomize cell collection before clusterization
    std::shuffle(cells.begin(), cells.end(), rnd);

    // clusterize assuming an 8-fold connectivity
    InPlaceClusterization::clusterize<SORT_AXIS, std::vector<Cell>,
                                      unsigned int,
                                      InPlaceClusterization::k8FoldConnection>(
        cells);

    // compute cluster hashes using the in-place clustered cell collection
    std::vector<std::size_t> cluster_hashes = computeClusterHashes(cells);

    // ... and test that input and output hashes agree.
    std::ranges::sort(cls, clHashComp);
    std::ranges::sort(cluster_hashes);

    BOOST_CHECK_EQUAL(cls.size(), cluster_hashes.size());
    for (std::size_t i = 0; i < cls.size(); i++) {
      BOOST_CHECK_EQUAL(cls.at(i).hash, cluster_hashes.at(i));
    }
  }
}

// Test for the in-place clusterization algorithm using sample cluster data
BOOST_AUTO_TEST_CASE(InPlaceClusterization_ClusterizeTestData) {
  // the sample cluster data
  std::vector<std::pair<bool, std::int16_t*>> test_data_sets{
      std::make_pair(false, cluster_test_data),
      std::make_pair(true, cluster_test_data8),
      std::make_pair(true, cluster_test_data8_2)};

  for (auto [fold8, test_data_input] : test_data_sets) {
    static constexpr unsigned int SORT_AXIS = 0;
    using Cell = InPlaceClusterization::Cell<std::int16_t, 2, std::uint16_t>;
    using Cluster = Cluster2D;

    // create a cluster collection from the sample cluster data
    std::vector<Cluster> cls = createClusterTestData(test_data_input);
    std::vector<Cell> cells;

    // create a cell collection
    unsigned int cell_idx = 0;
    for (const Cluster& cl : cls) {
      for (const Cell2D& a_cell : cl.cells) {
        // this test will not work if the coordinates or the index exceed
        // 16bit
        assert(static_cast<std::int16_t>(getCellRow(a_cell)) ==
               getCellRow(a_cell));
        assert(static_cast<std::int16_t>(getCellColumn(a_cell)) ==
               getCellColumn(a_cell));
        assert(static_cast<std::uint16_t>(cell_idx) == cell_idx);
        cells.push_back(Cell(
            std::array<std::int16_t, 2>{
                static_cast<std::int16_t>(getCellRow(a_cell)),
                static_cast<std::int16_t>(getCellColumn(a_cell))},
            static_cast<std::uint16_t>(cell_idx)));
        ++cell_idx;
      }
    }

    // shuffle the cell collection before testing the clusterization
    std::size_t startSeed = 71902647;
    std::mt19937_64 rnd(startSeed++);
    std::shuffle(cells.begin(), cells.end(), rnd);

    // clusterize using the in-place clusterization algorithm
    if (fold8) {
      InPlaceClusterization::clusterize<
          SORT_AXIS, std::vector<Cell>, unsigned int,
          InPlaceClusterization::k8FoldConnection>(cells);
    } else {
      InPlaceClusterization::clusterize<
          SORT_AXIS, std::vector<Cell>, unsigned int,
          InPlaceClusterization::k4FoldConnection>(cells);
    }

    // compute the per cluster hashes, using the in-place clustered cell
    // collection
    std::vector<std::size_t> cluster_hashes = computeClusterHashes(cells);

    // input cluster collection and the output cluster hashes
    std::ranges::sort(cls, clHashComp);
    std::ranges::sort(cluster_hashes);

    // and compare the input and output agree.
    BOOST_CHECK_EQUAL(cls.size(), cluster_hashes.size());
    for (std::size_t i = 0; i < cls.size(); i++) {
      BOOST_CHECK_EQUAL(cls.at(i).hash, cluster_hashes.at(i));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
