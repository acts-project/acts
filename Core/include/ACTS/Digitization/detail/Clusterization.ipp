// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

template <typename Cell>
std::vector<std::vector<Cell>>
Acts::createClusters(std::unordered_map<size_t, std::pair<Cell, bool>>& cellMap,
                     size_t nBins0,
                     bool   commonCorner,
                     double energyCut)
{
  // the output
  std::vector<std::vector<Cell>> mergedCells;
  // now go through cells and label
  for (auto& cell : cellMap) {
    // check if the cell was already used
    if (!(cell.second.second)
        && (cell.second.first.depositedEnergy() >= energyCut)) {
      // create new cluster
      mergedCells.push_back(std::vector<Cell>());
      // fill all cells belonging to that cluster
      ccl(mergedCells, cellMap, cell.first, nBins0, commonCorner, energyCut);
    }
  }
  // return the grouped together cells
  return mergedCells;
}

template <typename Cell>
void
Acts::ccl(std::vector<std::vector<Cell>>& mergedCells,
          std::unordered_map<size_t, std::pair<Cell, bool>>& cellMap,
          size_t index,
          size_t nBins0,
          bool   commonCorner,
          double energyCut)
{
  // add current cell to cluster
  auto cellA = cellMap.at(index).first;
  // check if cell energy is higher than energy threshold to activate the cell
  if (cellA.depositedEnergy() >= energyCut) {
    // add current cell to current cluster
    mergedCells.back().push_back(cellA);
    cellMap.at(index).second = true;
    // go recursively through all neighbours of this cell, if present
    // calculate neighbour indices first
    int              iMin = -1;
    int              jMin = -nBins0;
    int              iMax = 1;
    int              jMax = nBins0;
    std::vector<int> neighbourIndices;
    // the neighbour indices - filled depending on merging case
    if ((index % nBins0) == 0) {
      // left edge case
      if (commonCorner) {
        neighbourIndices = {jMin, jMin + iMax, iMax, jMax, jMax + iMax};
      } else {
        neighbourIndices = {jMin, iMax, jMax};
      }
    } else if (((index + 1) % nBins0) == 0) {
      // right edge case
      if (commonCorner) {
        neighbourIndices = {jMin + iMin, jMin, iMin, jMax + iMin, jMax};
      } else {
        neighbourIndices = {jMin, iMin, jMax};
      }
    } else {
      if (commonCorner) {
        neighbourIndices = {jMin + iMin,
                            jMin,
                            jMin + iMax,
                            iMin,
                            iMax,
                            jMax + iMin,
                            jMax,
                            jMax + iMax};
      } else {
        neighbourIndices = {jMin, iMin, iMax, jMax};
      }
    }
    // go through neighbours and recursively call connected component algorithm
    for (auto& i : neighbourIndices) {
      // calculate neighbour index of current cell
      int neighbourIndex = int(index) + i;
      // check if neighbour is there
      ///   auto startSearchMap = std::chrono::system_clock::now();
      auto search = cellMap.find(neighbourIndex);
      if ((search != cellMap.end())) {
        // get the corresponding index and call function again
        auto newIndex = search->first;
        if (!cellMap.at(newIndex).second) {
          ccl(mergedCells, cellMap, newIndex, nBins0, commonCorner, energyCut);
        }  // check if was used already
      }    // check if neighbour is there
    }      // go through neighbour indics
  }        // check energy cut
}