// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename external_spacepoint_t>
Acts::BinFinder<external_spacepoint_t>::BinFinder(
    const std::vector<std::pair<int, int>>& zBinNeighbors, int numPhiNeighbors)
    : m_zBinNeighbors(zBinNeighbors), m_numPhiNeighbors(numPhiNeighbors) {}

template <typename external_spacepoint_t>
boost::container::small_vector<size_t, 9>
Acts::BinFinder<external_spacepoint_t>::findBins(
    size_t phiBin, size_t zBin,
    const Acts::SpacePointGrid<external_spacepoint_t>* binnedSP) const {
  // if zBinNeighbors is not defined, get the indices using
  // neighborHoodIndices
  if (m_zBinNeighbors->empty()) {
    return binnedSP->neighborHoodIndices({phiBin, zBin}).collect();
  }
  // if the zBinNeighbors is defined, get the indices from there
  std::array<std::pair<int, int>, 2> sizePerAxis;
  sizePerAxis.at(0) = std::make_pair(-m_numPhiNeighbors, m_numPhiNeighbors);
  sizePerAxis.at(1) = (*m_zBinNeighbors)[zBin - 1];
  return binnedSP->neighborHoodIndices({phiBin, zBin}, sizePerAxis).collect();
}
