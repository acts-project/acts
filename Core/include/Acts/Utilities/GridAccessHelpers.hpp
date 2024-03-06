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
#include <tuple>
#include <vector>

namespace Acts {
namespace GridAccessHelpers {

/// @brief The local cast maps the local position on the surface where
/// the binning is defined to the actual lookup position of the grid.
///
/// Eg. for a cartesian 2 dimensional this would just take x and y, as they
/// map directly onto the grid structure. For a cylindrical surface with a z/phi
/// grid, however the local accesor needs to divide by the radius
///
/// @brief Shift, Scale and Access
///
/// @note this can be used to access from an r/phi local position
/// @note this can be used to shift the surface local frame to the grid frame
struct LocalAccess {
  std::size_t localIndex = 0u;
  ActsScalar shift = 0.;
  ActsScalar scale = 1.;

  /// @brief From local position to grid local axis
  /// @param lposition  the local position at input
  /// @return a scalar value in the grid local axis
  ActsScalar toGridLocal(const Acts::Vector2& lposition) const {
    return scale * (lposition[localIndex] + shift);
  }

  /// @brief From grid local to frame local
  /// @param v  the grid local axis value
  ActsScalar toFrameLocal(ActsScalar v) const { return v / scale - shift; }
};

/// Unroll the cast loop - from local position
///
/// @tparam cast_container is the container type of cast objects
/// @tparam Array is the array type to be filled
///
/// @param lposition is the local position at input
/// @param globalCasts is the cast value vector from global to grid position
/// @param ra is the array to be filled
template <typename local_cast_container, typename Array, std::size_t... idx>
void fillCasts(const Vector2& lposition, const local_cast_container& localCasts,
               Array& ra, std::index_sequence<idx...> /*indices*/) {
  ((ra[idx] = localCasts[idx].toGridLocal(lposition)), ...);
}

/// Cast into a lookup position - from local position
///
/// This method allows to transform a local position into a grid position
/// that is specified by LocalAccess structs.
///
/// @tparam local_cast_container is the container type of cast objects
/// @tparam Array is the array type to be filled
///
/// @param lposition is the position at input
/// @param localCasts is the cast value vector from global to grid position
///
/// @return a grid point in an appropriate format
template <typename grid_type, typename local_cast_container>
typename grid_type::point_t castLocal(const Vector2& lposition,
                                      const local_cast_container& localCasts) {
  // Fill the grid point from global
  typename grid_type::point_t casted{};
  fillCasts(lposition, localCasts, casted,
            std::make_integer_sequence<std::size_t, grid_type::DIM>{});
  return casted;
}

/// Unroll the cast loop - from global position
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

}  // namespace GridAccessHelpers
}  // namespace Acts
