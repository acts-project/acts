// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BinnedGroup.hpp"

#include <numeric>
#include <stdexcept>
#include <string>

namespace Acts {

namespace detail {

/// Validate a navigation array and fill in the axes that were left undefined.
///
/// Entries are 1-based local bin indices; 0 is the underflow bin, which is
/// never filled. An empty vector is replaced by 1, 2, ..., nBins.
///
/// @param navigation The navigation array, modified in place
/// @param numLocalBins Number of local bins per axis
/// @throw std::invalid_argument if an entry is out of range or repeated
template <std::size_t DIM>
void validateAndCompleteNavigationArray(
    std::array<std::vector<std::size_t>, DIM>& navigation,
    const std::array<std::size_t, DIM>& numLocalBins) {
  for (std::size_t i(0ul); i < DIM; ++i) {
    std::vector<std::size_t>& bins = navigation[i];
    const std::size_t nBins = numLocalBins[i];

    /// Undefined for this axis: default to a std::iota from 1ul
    if (bins.empty()) {
      bins.resize(nBins);
      std::iota(bins.begin(), bins.end(), 1ul);
      continue;
    }

    std::vector<bool> visited(nBins + 1ul, false);
    for (std::size_t bin : bins) {
      if (bin == 0ul || bin > nBins) {
        throw std::invalid_argument(
            "Invalid navigation for axis " + std::to_string(i) + ": bin " +
            std::to_string(bin) +
            " is out of range. Local bin indices are 1-based and must lie "
            "within [1, " +
            std::to_string(nBins) + "].");
      }
      if (visited[bin]) {
        throw std::invalid_argument(
            "Invalid navigation for axis " + std::to_string(i) + ": bin " +
            std::to_string(bin) + " is listed more than once.");
      }
      visited[bin] = true;
    }
  }
}

}  // namespace detail

template <typename grid_t>
BinnedGroup<grid_t>::BinnedGroup(
    grid_t&& grid, const GridBinFinder<BinnedGroup<grid_t>::DIM>& bottomFinder,
    const GridBinFinder<BinnedGroup<grid_t>::DIM>& topFinder,
    std::array<std::vector<std::size_t>, BinnedGroup<grid_t>::DIM> navigation)
    : m_grid(std::move(grid)),
      m_mask(m_grid.size(true), true),
      m_bottomBinFinder(&bottomFinder),
      m_topBinFinder(&topFinder),
      m_bins(std::move(navigation)) {
  detail::validateAndCompleteNavigationArray<DIM>(
      m_bins, m_grid.multiAxis().getNBins());
}

template <typename grid_t>
BinnedGroup<grid_t>::BinnedGroup(
    grid_t&& grid, std::vector<bool> mask,
    const GridBinFinder<BinnedGroup<grid_t>::DIM>& bottomFinder,
    const GridBinFinder<BinnedGroup<grid_t>::DIM>& topFinder,
    std::array<std::vector<std::size_t>, BinnedGroup<grid_t>::DIM> navigation)
    : m_grid(std::move(grid)),
      m_mask(std::move(mask)),
      m_bottomBinFinder(&bottomFinder),
      m_topBinFinder(&topFinder),
      m_bins(std::move(navigation)) {
  // Check the elements in the mask corresponds to all the global bins in the
  // grid so that we can check if a global bin is masked
  if (m_mask.size() != m_grid.size(true)) {
    throw std::invalid_argument(
        "Provided mask does not match the grid. The number of entries must "
        "correspond to the number of global bins in the grid.");
  }

  detail::validateAndCompleteNavigationArray<DIM>(
      m_bins, m_grid.multiAxis().getNBins());
}

template <typename grid_t>
const grid_t& BinnedGroup<grid_t>::grid() const {
  return m_grid;
}

template <typename grid_t>
grid_t& BinnedGroup<grid_t>::grid() {
  return m_grid;
}

template <typename grid_t>
const std::vector<bool>& BinnedGroup<grid_t>::mask() const {
  return m_mask;
}

template <typename grid_t>
BinnedGroupIterator<grid_t> BinnedGroup<grid_t>::begin() const {
  return BinnedGroupIterator<grid_t>(
      *this, std::array<std::size_t, BinnedGroup<grid_t>::DIM>(), m_bins);
}

template <typename grid_t>
BinnedGroupIterator<grid_t> BinnedGroup<grid_t>::end() const {
  std::array<std::size_t, BinnedGroup<grid_t>::DIM> endline{};
  for (std::size_t i(0ul); i < BinnedGroup<grid_t>::DIM; ++i) {
    endline[i] = m_bins[i].size();
  }
  return BinnedGroupIterator<grid_t>(*this, endline, m_bins);
}

}  // namespace Acts
