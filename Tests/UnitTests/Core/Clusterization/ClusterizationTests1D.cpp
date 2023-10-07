// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Clusterization/Clusterization.hpp"

#include <boost/functional/hash.hpp>

namespace Acts {
namespace Test {

struct Cell1D {
  Cell1D(int colv) : col(colv) {}
  int col;
  Ccl::Label label{Ccl::NO_LABEL};
};

bool cellComp(const Cell1D& left, const Cell1D& right) {
  return left.col < right.col;
}

Ccl::Label& getCellLabel(Cell1D& cell) {
  return cell.label;
}

int getCellColumn(const Cell1D& cell) {
  return cell.col;
}

struct Cluster1D {
  std::vector<Cell1D> cells;
  size_t hash{};
};

void clusterAddCell(Cluster1D& cl, const Cell1D& cell) {
  cl.cells.push_back(cell);
}

bool clHashComp(const Cluster1D& left, const Cluster1D& right) {
  return left.hash < right.hash;
}

void hash(Cluster1D& cl) {
  std::sort(cl.cells.begin(), cl.cells.end(), cellComp);
  cl.hash = 0;
  for (const Cell1D& c : cl.cells) {
    boost::hash_combine(cl.hash, c.col);
  }
}

BOOST_AUTO_TEST_CASE(Grid_1D_rand) {
  using Cell = Cell1D;
  using CellC = std::vector<Cell>;
  using Cluster = Cluster1D;
  using ClusterC = std::vector<Cluster>;

  size_t minsize = 1;
  size_t maxsize = 10;
  size_t minspace = 1;
  size_t maxspace = 10;
  size_t nclusters = 100;
  int startSeed = 204769;
  size_t ntries = 100;

  std::cout << "Grid_1D_rand test with parameters:" << std::endl;
  std::cout << "  minsize = " << minsize << std::endl;
  std::cout << "  maxsize = " << maxsize << std::endl;
  std::cout << "  minspace = " << minspace << std::endl;
  std::cout << "  maxspace = " << maxspace << std::endl;
  std::cout << "  nclusters = " << nclusters << std::endl;
  std::cout << "  startSeed = " << startSeed << std::endl;
  std::cout << "  ntries = " << ntries << std::endl;

  while (ntries-- > 0) {
    std::mt19937_64 rnd(startSeed++);
    std::uniform_int_distribution<size_t> distr_size(minsize, maxsize);
    std::uniform_int_distribution<size_t> distr_space(minspace, maxspace);

    int col = 0;

    CellC cells;
    ClusterC clusters;
    for (size_t i = 0; i < nclusters; i++) {
      Cluster cl;
      col += distr_space(rnd);
      size_t size = distr_size(rnd);
      for (size_t j = 0; j < size; j++) {
        Cell cell(col++);
        cells.push_back(cell);
        clusterAddCell(cl, cell);
      }
      clusters.push_back(std::move(cl));
    }
    for (Cluster& cl : clusters) {
      hash(cl);
    }

    std::shuffle(cells.begin(), cells.end(), rnd);

    ClusterC newCls = Ccl::createClusters<CellC, ClusterC, 1>(cells);

    for (Cluster& cl : newCls) {
      hash(cl);
    }

    std::sort(clusters.begin(), clusters.end(), clHashComp);
    std::sort(newCls.begin(), newCls.end(), clHashComp);

    BOOST_CHECK_EQUAL(clusters.size(), newCls.size());
    for (size_t i = 0; i < clusters.size(); i++) {
      BOOST_CHECK_EQUAL(clusters.at(i).hash, newCls.at(i).hash);
    }
  }
}

}  // namespace Test
}  // namespace Acts
