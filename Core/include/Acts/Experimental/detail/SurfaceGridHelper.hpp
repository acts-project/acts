// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Enumerate.hpp"

#include <algorithm>
#include <vector>

namespace Acts {

namespace Experimental {
namespace detail {

/// @brief Checks if global bin is valid
/// @param grid is the grid for the check
/// @param bin the global bin index
/// @return bool if the bin is valid
/// @note Valid means that the index points to a bin which is not a under
///       or overflow bin or out of range in any axis.
template <typename Grid_t>
static bool isValidBin(Grid_t& grid, std::size_t bin) {
  const auto axes = grid.axes();
  std::array<std::size_t, axes.size()> indices =
      grid.localBinsFromGlobalBin(bin);
  std::array<std::size_t, axes.size()> nBins = grid.numLocalBins();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    auto idx = indices.at(i);
    if (idx <= 0 || idx >= nBins.at(i) + 1) {
      return false;
    }
  }
  return true;
}

template <typename Grid_t>
static void populateNeighborhood(Grid_t& grid) {
  // The cache
  std::map<size_t, typename Grid_t::value_type> neighborCache;
  // First loop w/o filling:
  // - calculate neighbors for every bin and store in map
  for (std::size_t i = 0u; i < grid.size(); i++) {
    if (isValidBin(grid, i)) {
      typename Grid_t::index_t loc = grid.localBinsFromGlobalBin(i);
      auto neighborIdxs = grid.neighborHoodIndices(loc, 1u);
      std::vector<std::size_t> neighbors;
      // Loop over the neightborhood indices
      for (const auto idx : neighborIdxs) {
        if (idx < grid.size()) {
          const auto& binContent = grid.at(idx);
          std::copy(binContent.begin(), binContent.end(),
                    std::back_inserter(neighbors));
        }
      }
      neighborCache[i] = neighbors;
    }
  }
  // Second loop w filling
  for (auto [i, nc] : neighborCache) {
    auto& binContent = grid.at(i);
    for (const auto& in : nc) {
      if (std::find(binContent.begin(), binContent.end(), in) ==
          binContent.end()) {
        binContent.push_back(in);
      }
    }
  }
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
