// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename external_spacepoint_t>
Acts::BinFinder<external_spacepoint_t>::BinFinder(
    const std::vector<std::vector<int> >&& indicesVector)
    : m_indicesVector(std::move(indicesVector)) {}

template <typename external_spacepoint_t>
boost::container::small_vector<size_t, 10>
Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) {
  boost::container::small_vector<size_t, 9> indices;
  // if the map is not defined, get the indices using neighborHoodIndices
  if (m_indicesVector.empty()) {
    indices = binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
  }
  // if the map is defined, get the indices from the map
  else {
    for (int phiBinIndex = -1; phiBinIndex <= +1; phiBinIndex++) {
      for (std::vector<int>::size_type zBinIndex = 0;
           zBinIndex < m_indicesVector[zBin - 1].size(); zBinIndex++) {
        unsigned long zBinGlobalIndex = m_indicesVector[zBin - 1][zBinIndex];
        auto phiBinGlobalIndex = phiBin + phiBinIndex;
        if (phiBinGlobalIndex == 0) {
          phiBinGlobalIndex = (binnedSP->numLocalBins())[0];
        } else if (phiBinGlobalIndex == 1 + (binnedSP->numLocalBins())[0]) {
          phiBinGlobalIndex = 1;
        }
        const std::array<size_t, 2> globalIndicesArray = {phiBinGlobalIndex,
                                                          zBinGlobalIndex};
        auto globalIndice =
            binnedSP->globalBinFromLocalBins(globalIndicesArray);
        indices.push_back(globalIndice);
      }
    }
  }

  return {indices.begin(), indices.end()};
}
