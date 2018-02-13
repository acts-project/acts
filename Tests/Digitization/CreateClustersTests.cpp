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
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/test/data/test_case.hpp>
#include <map>
#include <utility>
#include <vector>
#include "ACTS/Digitization/DigitizationCell.hpp"
#include "ACTS/Digitization/Segmentation.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
      Graph;

  const std::map<size_t, std::vector<Acts::DigitizationCell>>
  createClusters(const std::vector<Acts::DigitizationCell>& cells,
                 bool commonCorner = false)
  {
    // the output
    std::map<size_t, std::vector<Acts::DigitizationCell>> mergedCells;
    // create the graph
    Graph graph;
    //  add the needed amount of vertices to the graph
    for (auto cell : cells) {
      // add vertex to graph for each cell
      add_vertex(graph);
    }

    // loop through cells and set the edges
    for (auto cellA = std::begin(cells); cellA != (std::end(cells) - 1);
         cellA++) {
      for (auto cellB = (cellA + 1); cellB != std::end(cells); cellB++) {
        // the channels
        auto c0A = cellA->channel0;
        auto c1A = cellA->channel1;
        auto c0B = cellB->channel0;
        auto c1B = cellB->channel1;
        // the conditions
        bool isNeighbour0   = (c0A == (c0B - 1) || c0A == (c0B + 1));
        bool isSameChannel0 = (c0A == c0B);
        bool isNeighbour1   = (c1A == (c1B - 1) || c1A == (c1B + 1));
        bool isSameChannel1 = (c1A == c1B);
        // distinguish between between merging cells, when they have a common
        // corner (8cell) or a common edge (4cell)
        if (commonCorner) {
          // 8cell
          if ((isNeighbour0 || isSameChannel0)
              && (isNeighbour1 || isSameChannel1)) {
            auto indexA = std::distance(cells.begin(), cellA);
            add_edge(indexA, indexA + 1, graph);
          }
        } else {
          // 4cell
          if ((isNeighbour0 && isSameChannel1)
              || (isNeighbour1 && isSameChannel0)) {
            auto indexA = std::distance(cells.begin(), cellA);
            add_edge(indexA, indexA + 1, graph);
          }
        }
      }  // cellsB
    }    // cellsA

    std::vector<size_t> component(num_vertices(graph));
    size_t              num = connected_components(graph, &component[0]);

    std::cout << "Total number of components: " << num << std::endl;
    for (size_t i = 0; i != component.size(); ++i) {
      mergedCells[component[i]].push_back(cells.at(i));
      std::cout << "Vertex " << i << " is in component " << component[i]
                << std::endl;
    }
    return mergedCells;
  }

  BOOST_AUTO_TEST_CASE(merge_clusters)
  {
    std::vector<Acts::DigitizationCell> testCells1;
    testCells1.push_back(Acts::DigitizationCell(2, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(2, 4, 1));
    testCells1.push_back(Acts::DigitizationCell(3, 3, 1));
    testCells1.push_back(Acts::DigitizationCell(8, 9, 1));

    auto mergedCells1 = createClusters(testCells1);
    std::cout << "#clusters: " << mergedCells1.size() << std::endl;
    for (auto& a : mergedCells1) {
      std::cout << "# of cells in cluster " << a.first << ": "
                << a.second.size() << std::endl;
    }
    std::cout << "New Test" << std::endl;
    auto mergedCells2 = createClusters(testCells1, true);
    std::cout << "#clusters: " << mergedCells2.size() << std::endl;
    for (auto& a : mergedCells2) {
      std::cout << "# of cells in cluster " << a.first << ": "
                << a.second.size() << std::endl;
    }
  }
}
}
