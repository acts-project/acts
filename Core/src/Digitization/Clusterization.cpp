// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Digitization/Clusterization.hpp"
#include <iostream>

std::vector<Acts::DigitizationCell>
Acts::mergeCells(std::vector<Acts::DigitizationCell>& cells,
                 bool                                 analogueReadout,
                 double                               energyCut)
{
  if (!cells.size()) return cells;
  // the output
  std::vector<std::vector<Acts::DigitizationCell>> mergedCells;
  // create the graph
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
  //  add the needed amount of vertices to the graph
  for (auto& cell : cells) {
    (void)cell;
    // add vertex to graph for each cell
    add_vertex(graph);
  }

  // loop through cells and set the edges
  for (auto cellA = std::begin(cells); cellA != (std::end(cells) - 1);
       cellA++) {
    auto indexA = std::distance(cells.begin(), cellA);
    for (auto cellB = (cellA + 1); cellB != std::end(cells); cellB++) {
      auto indexB = std::distance(cells.begin(), cellB);
      // the channels
      auto c0A = cellA->channel0;
      auto c1A = cellA->channel1;
      auto c0B = cellB->channel0;
      auto c1B = cellB->channel1;

      // the conditions
      bool isSameChannel0 = (c0A == c0B);
      bool isSameChannel1 = (c1A == c1B);
      // same cell
      if (isSameChannel0 && isSameChannel1) {
        add_edge(indexA, indexB, graph);
      }  // if
    }    // cellsB
  }      // cellsA

  std::vector<size_t> component(num_vertices(graph));
  connected_components(graph, &component[0]);
  // copy the component map
  std::vector<size_t> keys = component;
  // sort
  std::sort(keys.begin(), keys.end());
  // get the keys
  auto it = std::unique(keys.begin(), keys.end());
  keys.erase(it, keys.end());

  for (auto& comp : keys) {
    std::vector<Acts::DigitizationCell> compCells;
    for (size_t i = 0; i != component.size(); ++i) {
      if (component[i] == comp) {
        compCells.push_back(cells.at(i));
      }
    }
    mergedCells.push_back(compCells);
  }
  std::vector<Acts::DigitizationCell> sameCells;

  for (auto& sCells : mergedCells) {
    Acts::DigitizationCell& newCell = *std::begin(sCells);
    for (auto sCell = (std::begin(sCells) + 1); sCell != std::end(sCells);
         sCell++) {
      if (analogueReadout)
        newCell.data += sCell->depositedEnergy(analogueReadout);
    }
    // only add cell if is above threshold
    if (newCell.depositedEnergy(analogueReadout) >= energyCut)
      sameCells.push_back(newCell);
  }

  return sameCells;
}

std::vector<std::vector<Acts::DigitizationCell>>
Acts::createClusters(const std::vector<Acts::DigitizationCell>& cells,
                     bool                                       commonCorner)
{
  // the output
  std::vector<std::vector<Acts::DigitizationCell>> mergedCells;
  if (!cells.size()) return mergedCells;

  // create the graph
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
  //  add the needed amount of vertices to the graph
  for (auto& cell : cells) {
    (void)cell;
    // add vertex to graph for each cell
    add_vertex(graph);
  }

  // loop through cells and set the edges
  for (auto cellA = std::begin(cells); cellA != (std::end(cells) - 1);
       cellA++) {
    auto indexA = std::distance(cells.begin(), cellA);
    for (auto cellB = (cellA + 1); cellB != std::end(cells); cellB++) {
      auto indexB = std::distance(cells.begin(), cellB);
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
      // if the cell is at the same position they need to be merged
      // distinguish between between merging cells, when they have a common
      // corner (8cell) or a common edge (4cell)
      if (commonCorner) {
        // 8cell
        if ((isNeighbour0 || isSameChannel0)
            && (isNeighbour1 || isSameChannel1)) {
          add_edge(indexA, indexB, graph);
        }
      } else {
        // 4cell
        if ((isNeighbour0 && isSameChannel1)
            || (isNeighbour1 && isSameChannel0)) {
          add_edge(indexA, indexB, graph);
        }
      }
    }  // cellsB
  }    // cellsA

  std::vector<size_t> component(num_vertices(graph));
  connected_components(graph, &component[0]);
  // copy the component map
  std::vector<size_t> keys = component;
  // sort
  std::sort(keys.begin(), keys.end());
  // get the keys
  auto it = std::unique(keys.begin(), keys.end());
  keys.erase(it, keys.end());

  for (auto& comp : keys) {
    std::vector<Acts::DigitizationCell> compCells;
    for (size_t i = 0; i != component.size(); ++i) {
      if (component[i] == comp) {
        compCells.push_back(cells.at(i));
      }
    }
    mergedCells.push_back(compCells);
  }
  return mergedCells;
}
