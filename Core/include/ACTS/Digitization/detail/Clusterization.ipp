// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>

template <typename Cell>
std::vector<std::vector<Cell>>
Acts::createClusters(const std::vector<Cell>& cells,
                     size_t                   nBins0,
                     size_t                   nBins1,
                     bool                     commonCorner,
                     bool                     analogueReadout,
                     double                   energyCut)
{
  // the output
  std::vector<std::vector<Cell>> mergedCells;
  // map containing all activated cells + boolean indicating if they have been
  // used already
  std::unordered_map<size_t, std::pair<Cell, bool>> cellMap;
  // Now fill cell map and merge cells if requested
  for (auto& cell : cells) {
    // calculate key of cell which is the global grid index
    size_t globalIndex = cell.channel0 + nBins0 * cell.channel1;
    // insert new cell - only adds cell, if cell was not there yet
    auto insertCell = cellMap.insert({globalIndex, {cell, false}});
    if (!insertCell.second) {
      // check if there is already a cell at same position and merge in that
      // case
      insertCell.first->second.first.addCell(cell, analogueReadout);
    }
  }
  // now go through cells and label
  for (auto& cell : cellMap) {
    // check if the cell was already used
    if (!(cell.second.second)
        && (cell.second.first.depositedEnergy(analogueReadout) >= energyCut)) {
      // create new cluster
      mergedCells.push_back(std::vector<Cell>());
      // fill all cells belonging to that cluster
      ccl(mergedCells,
          cellMap,
          cell.first,
          nBins0,
          nBins1,
          commonCorner,
          analogueReadout,
          energyCut);
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
          size_t nBins1,
          bool   commonCorner,
          bool   analogueReadout,
          double energyCut)
{
  // add current cell to cluster
  auto cellA = cellMap.at(index).first;
  // check if cell energy is higher than energy threshold to activate the cell
  if (cellA.depositedEnergy(analogueReadout) >= energyCut) {
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
    } else if (index == 0) {
      // left corner case
      if (commonCorner) {
        neighbourIndices = {iMax, jMax, jMax + iMax};
      } else {
        neighbourIndices = {iMax, jMax};
      }
    } else if (((index + 1) % nBins0) == 0) {
      // right edge case
      if (commonCorner) {
        neighbourIndices = {jMin + iMin, jMin, iMin, jMax + iMin, jMax};
      } else {
        neighbourIndices = {jMin, iMin, jMax};
      }
    } else if (index == (nBins0 * nBins1)) {
      // right corner case
      if (commonCorner) {
        neighbourIndices = {jMin + iMin, jMin, iMin};
      } else {
        neighbourIndices = {jMin, iMin};
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
      // check of neighbour is there
      auto search = cellMap.find(neighbourIndex);
      if ((search != cellMap.end())) {
        // get the corresponding index and call function again
        auto newIndex = search->first;
        if (!cellMap.at(newIndex).second) {
          ccl(mergedCells,
              cellMap,
              newIndex,
              nBins0,
              nBins1,
              commonCorner,
              analogueReadout,
              energyCut);
        }  // check if was used already
      }    // check if neighbour is there
    }      // go through neighbour indics
  }        // check energy cut
}
