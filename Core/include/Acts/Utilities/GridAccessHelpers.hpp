// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <stdexcept>
#include <utility>

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
/// @note void function that fills the provided array
template <typename Array, typename local_indices, std::size_t... idx>
void fillLocal(const Vector2& lposition, const local_indices& laccess,
               Array& ra, std::index_sequence<idx...> /*indices*/) {
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
template <typename grid_type, typename local_indices>
typename grid_type::point_t accessLocal(const Vector2& lposition,
                                        const local_indices& laccess) {
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

}  // namespace GridAccessHelpers
}  // namespace Acts
