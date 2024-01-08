// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <vector>

namespace Acts {
namespace GridAccessHelpers {

/// Unroll the cast loop
///
/// @tparam cast_container is the container type of cast objects
/// @tparam Array is the array type to be filled
///
/// @param position is the position of the update call
/// @param globalCasts is the cast value vector from global to grid position
/// @param ra is the array to be filled
template <typename cast_container, typename Array, std::size_t... idx>
void fillCasts(const Vector3& position, const cast_container& globalCasts,
               Array& ra, std::index_sequence<idx...> /*indices*/) {
  ((ra[idx] = VectorHelpers::cast(position, globalCasts[idx])), ...);
}
/// Cast into a lookup position
///
/// This method allows to transform a global position into a grid position
/// that is specified by binning values. It casts the global position info
/// into a format suitable for the grid.
///
/// @tparam cast_container is the container type of cast objects
/// @tparam Array is the array type to be filled
///
/// @param position is the position of the update call
/// @param globalCasts is the cast value vector from global to grid position
///
/// @return a grid point in an appropriate format
template <typename grid_type, typename cast_container>
typename grid_type::point_t castPosition(const Vector3& position,
                                         const cast_container& globalCasts) {
  // Fill the grid point from global
  typename grid_type::point_t casted{};
  fillCasts(position, globalCasts, casted,
            std::make_integer_sequence<std::size_t, grid_type::DIM>{});
  return casted;
}

/// Unroll the local position loop
///
/// @param lposition is the local position
/// @param laccess the local accessors
/// @param ra is the array to be filled
///
/// @note voide function that fill the provided array
template <typename Array, std::size_t... idx>
void fillLocal(const Vector2& lposition,
               const std::vector<std::size_t>& laccess, Array& ra,
               std::index_sequence<idx...> /*indices*/) {
  ((ra[idx] = lposition[laccess[idx]]), ...);
}

/// Access local parameters for a propriate lookup position
///
/// This method allows to transform a local position into a
/// grid position, it only works for 1-D and 2-D grids.
///
/// @param lposition is the position of the update call
/// @param laccess the local accessors
///
/// @return an array suitable for the grid
template <typename grid_type>
typename grid_type::point_t accessLocal(
    const Vector2& lposition, const std::vector<std::size_t>& laccess) {
  if constexpr (grid_type::DIM > 2u) {
    throw std::invalid_argument(
        "GridAccessHelper: only 1-D and 2-D grids are possible for local "
        "access.");
  }
  // Fill the grid point from local according to the accessors
  typename grid_type::point_t accessed{};
  fillLocal(lposition, laccess, accessed,
            std::make_integer_sequence<std::size_t, grid_type::DIM>{});
  return accessed;
}

/// Helper method to convert:
/// 1-D grids:
/// * [ 0 - n ]
/// * [ 0 - m ]
/// 2-D grid:
/// * [ 0 - n ] x [ 0, m ]
///
/// To a local position that can be looked up the accessLocal(...) method
///
/// @param grid is the grid, we need the axes
/// @param bin0 is the first bin
/// @param bin1 is the second bin
/// @param laccess the local accessors for 1D grid case
///
/// @note this method is not very fast, it shall not be used for performance
///      critical code
///
/// @return a local position that can be looked up via accessLocal(...)
template <typename grid_type>
Vector2 toLocal(const grid_type& grid, std::size_t bin0, std::size_t bin1,
                std::size_t laccess = 0u) {
  Vector2 lposition{};
  auto gridAxes = grid.axes();
  // One-dimensional case, needs decision which one is assigned
  if constexpr (grid_type::DIM == 1u) {
    if (laccess > 1u) {
      throw std::invalid_argument(
          "GridAccessHelper: only 0u/1u are allowed for local access.");
    }

    // Get axis for bin edges
    std::size_t bin = laccess == 0u ? bin0 : bin1;
    const auto& edges = gridAxes[laccess]->getBinEdges();
    ActsScalar pval = 0.5 * (edges[bin] + edges[bin + 1]);
    lposition[laccess] = pval;
  }
  // Two-dimensional case, relatively straight forward
  if constexpr (grid_type::DIM == 2u) {
    // Get axis for the bin edge
    const auto& edges0 = gridAxes[0u]->getBinEdges();
    const auto& edges1 = gridAxes[1u]->getBinEdges();
    ActsScalar pval0 = 0.5 * (edges0[bin0] + edges0[bin0 + 1u]);
    ActsScalar pval1 = 0.5 * (edges1[bin1] + edges1[bin1 + 1u]);
    lposition = {pval0, pval1};
  }

  // Error for DIM > 2
  if constexpr (grid_type::DIM > 2u) {
    throw std::invalid_argument(
        "GridAccessHelper: only 1-D and 2-D grids are possible for binned "
        "access.");
  }

  // Return the suitable local position
  return lposition;
}

}  // namespace GridAccessHelpers
}  // namespace Acts
