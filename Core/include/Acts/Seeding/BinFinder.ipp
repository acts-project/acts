// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename external_spacepoint_t>
Acts::BinFinder<external_spacepoint_t>::BinFinder(
    const std::vector<std::vector<size_t> >&& indicesVector,
    const size_t&& numPhiNeighbors)
    : m_indicesVector(std::move(indicesVector)),
      m_numPhiNeighbors(std::move(numPhiNeighbors)) {}

template <typename external_spacepoint_t>
boost::container::small_vector<size_t, 10>
Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) {
  boost::container::small_vector<size_t, 9> indices;
  // if indicesVector is not defined, get the indices using neighborHoodIndices
  if (m_indicesVector.empty()) {
    indices = binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
  }
  // if the indicesVector is defined, get the indices from there
  else {
    // loop over the phi range defined by m_numPhiNeighbors
    int phiNeighborRange = (m_numPhiNeighbors - 1) / 2;
    for (int phiBinIndex = -phiNeighborRange; phiBinIndex <= phiNeighborRange;
         phiBinIndex++) {
      // loop over the z bins inside indicesVector
      for (size_t zBinIndex = 0; zBinIndex < m_indicesVector[zBin - 1].size();
           zBinIndex++) {
        // get z bin local index from indicesVector
        auto zBinLocalIndex = m_indicesVector[zBin - 1][zBinIndex];
        // get phi bin local index
        int phiBinLocalIndex = phiBin + phiBinIndex;
        int maxPhiBin = (binnedSP->numLocalBins())[0];
        // check if phiBin + phiBinIndex is outside of the grid: If phiBin +
        // phiBinIndex is smalle than the first bin, phiBinLocalIndex goes to
        // the end of the grid. If phiBin + phiBinIndex is greater than the last
        // bin, phiBinLocalIndex goes back to the first bins of the grid
        if (phiBinLocalIndex < 1) {
          phiBinLocalIndex += maxPhiBin;
        } else if (phiBinLocalIndex > maxPhiBin) {
          phiBinLocalIndex -= maxPhiBin;
        }
        const std::array<size_t, 2> localIndexArray = {phiBinLocalIndex,
                                                       zBinLocalIndex};
        // get the global bin index from local bin index
        auto globalIndex = binnedSP->globalBinFromLocalBins(localIndexArray);
        indices.push_back(globalIndex);
      }
    }
  }

  return {indices.begin(), indices.end()};
}
