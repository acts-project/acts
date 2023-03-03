// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <iosfwd>
#include <limits>

namespace Acts {

/// Tolerance for being numerical equal for geometry building
static constexpr ActsScalar s_epsilon =
    3 * std::numeric_limits<ActsScalar>::epsilon();

/// Tolerance for being on Surface
///
/// @note This is intentionally given w/o an explicit unit to avoid having
///       to include the units header unneccessarily. With the native length
///       unit of mm this corresponds to 0.1um.
static constexpr ActsScalar s_onSurfaceTolerance = 1e-4;

/// Tolerance for not being within curvilinear projection
/// this allows using the same curvilinear frame to eta = 6,
/// validity tested with IntegrationTests/PropagationTest
static constexpr ActsScalar s_curvilinearProjTolerance = 0.999995;

/// @enum NavigationDirection
/// The navigation direction is always with
/// respect to a given momentum or direction
enum class NavigationDirection : int { Backward = -1, Forward = 1 };

/// Convert navigation dir to index [0,1] which allows to
/// store direction dependent objects in std::array<T,2u>
///
/// @param nDir is the navigation direction at input
///
/// returns either 0 or 1
inline constexpr size_t indexFromDirection(NavigationDirection nDir) {
  if (nDir == NavigationDirection::Backward) {
    return 0u;
  }
  return 1u;
}

/// Convert and ndex [0,1] to a navigation direction
/// for sorting  in std::array<T,2u>
///
/// @param index is the navigation direction at input
///
/// returns either 0 or 1
inline constexpr NavigationDirection directionFromIndex(size_t index) {
  if (index == 0u) {
    return NavigationDirection::Backward;
  }
  return NavigationDirection::Forward;
}

/// This turns a signed value into a navigation direction
///
/// @param value is the signed value
///
/// @return a navigation direciton enum
inline constexpr NavigationDirection directionFromStepSize(ActsScalar value) {
  assert(value != 0);
  return value > 0 ? NavigationDirection::Forward
                   : NavigationDirection::Backward;
}

/// Invert a navigation direction enum
///
/// @param nDir is the navigation direction at input
///
/// return an opposite navigation direction
inline constexpr NavigationDirection invertDirection(NavigationDirection nDir) {
  return (nDir == NavigationDirection::Forward) ? NavigationDirection::Backward
                                                : NavigationDirection::Forward;
}

std::ostream& operator<<(std::ostream& os, NavigationDirection navDir);

// NavigationDirection * T

inline constexpr auto operator*(NavigationDirection dir, int value) {
  return static_cast<std::underlying_type_t<NavigationDirection>>(dir) * value;
}

inline constexpr auto operator*(NavigationDirection dir, float value) {
  return static_cast<std::underlying_type_t<NavigationDirection>>(dir) * value;
}

inline constexpr auto operator*(NavigationDirection dir, double value) {
  return static_cast<std::underlying_type_t<NavigationDirection>>(dir) * value;
}

inline Acts::Vector3 operator*(NavigationDirection dir,
                               const Acts::Vector3& value) {
  return static_cast<std::underlying_type_t<NavigationDirection>>(dir) * value;
}

// T * NavigationDirection

inline constexpr auto operator*(int value, NavigationDirection dir) {
  return value * static_cast<std::underlying_type_t<NavigationDirection>>(dir);
}

inline constexpr auto operator*(float value, NavigationDirection dir) {
  return value * static_cast<std::underlying_type_t<NavigationDirection>>(dir);
}

inline constexpr auto operator*(double value, NavigationDirection dir) {
  return value * static_cast<std::underlying_type_t<NavigationDirection>>(dir);
}

inline Acts::Vector3 operator*(const Acts::Vector3& value,
                               NavigationDirection dir) {
  return value * static_cast<std::underlying_type_t<NavigationDirection>>(dir);
}

// T *= NavigationDirection

inline constexpr auto operator*=(int& value, NavigationDirection dir) {
  value *= static_cast<std::underlying_type_t<NavigationDirection>>(dir);
  return value;
}

inline constexpr auto operator*=(float& value, NavigationDirection dir) {
  value *= static_cast<std::underlying_type_t<NavigationDirection>>(dir);
  return value;
}

inline constexpr auto operator*=(double& value, NavigationDirection dir) {
  value *= static_cast<std::underlying_type_t<NavigationDirection>>(dir);
  return value;
}

inline Acts::Vector3& operator*=(Acts::Vector3& value,
                                 NavigationDirection dir) {
  value *= static_cast<std::underlying_type_t<NavigationDirection>>(dir);
  return value;
}

///  This is a steering enum to tell which material update stage:
/// - PreUpdate  : update on approach of a surface
/// - FullUpdate : update when passing a surface
/// - PostUpdate : update when leaving a surface
enum class MaterialUpdateStage : int {
  PreUpdate = -1,
  FullUpdate = 0,
  PostUpdate = 1
};

std::ostream& operator<<(std::ostream& os, MaterialUpdateStage matUpdate);

/// @enum NoiseUpdateMode to tell how to deal with noise term in covariance
/// transport
/// - removeNoise: subtract noise term
/// - addNoise: add noise term
enum NoiseUpdateMode : int { removeNoise = -1, addNoise = 1 };

/// Components of coordinate vectors.
///
/// To be used to access coordinate components by named indices instead of magic
/// numbers. This must be a regular `enum` and not a scoped `enum class` to
/// allow implicit conversion to an integer. The enum value are thus visible
/// directly in `namespace Acts`.
///
/// This index enum is not user-configurable (in contrast e.g. to the track
/// parameter index enums) since it must be compatible with varying
/// dimensionality (2d-4d) and other access methods (`.{x,y,z}()` accessors).
enum CoordinateIndices : unsigned int {
  // generic position-like access
  ePos0 = 0,
  ePos1 = 1,
  ePos2 = 2,
  eTime = 3,
  // generic momentum-like access
  eMom0 = ePos0,
  eMom1 = ePos1,
  eMom2 = ePos2,
  eEnergy = eTime,
  // Cartesian spatial coordinates
  eX = ePos0,
  eY = ePos1,
  eZ = ePos2,
};

}  // namespace Acts
