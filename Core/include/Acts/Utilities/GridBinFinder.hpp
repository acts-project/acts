// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"

#include <variant>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @class BinFinder
/// The BinFinder is used by the ISPGroupSelector. It can be
/// used to find both bins that could be bottom bins as well as bins that could
/// be top bins, which are assumed to be the same bins. Does not take
/// interaction region into account to limit z-bins.
template <std::size_t DIM>
class GridBinFinder {
 public:
  template <typename... args>
  GridBinFinder(args&&... vals);

  /// Return all bins that could contain space points that can be used with the
  /// space points in the bin with the provided indices to create seeds.
  /// @param phiBin phi index of bin with middle space points
  /// @param zBin z index of bin with middle space points
  /// @param binnedSP phi-z grid containing all bins
  template <typename stored_t, class... Axes>
  boost::container::small_vector<std::size_t, Acts::detail::ipow(3, DIM)> findBins(
      const std::array<std::size_t, DIM>& locPosition,
      const Acts::Grid<stored_t, Axes...>& grid) const;

 private:
  template <typename first_value_t, typename... vals>
  void storeValue(first_value_t&& fv, vals&&... others);

  std::array<std::pair<int, int>, DIM> getSizePerAxis(
      const std::array<std::size_t, DIM>& locPosition) const;

 private:
  using stored_values_t = std::variant<int, std::vector<std::pair<int, int>>>;
  std::array<stored_values_t, DIM> m_values{};
};

}  // namespace Acts
#include "Acts/Utilities/GridBinFinder.ipp"
