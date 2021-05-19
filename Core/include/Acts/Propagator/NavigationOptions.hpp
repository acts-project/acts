// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Units.hpp"

#include <limits>

namespace Acts {

/// @brief struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The navigation direction
  NavigationDirection navDir = forward;

  /// The boundary check directive
  BoundaryCheck boundaryCheck = true;

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// object to check against: at start
  const object_t* startObject = nullptr;
  /// object to check against: at end
  const object_t* endObject = nullptr;

  /// Target surface to exclude
  const Surface* targetSurface = nullptr;
  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The maximum path limit for this navigation step
  double pathLimit = std::numeric_limits<double>::max();

  /// The overstep tolerance for this navigation step
  /// @note must be negative as it describes overstepping
  /// @todo could be dynamic in the future (pT dependent)
  double overstepLimit = -1 * UnitConstants::um;

  /// Constructor
  ///
  /// @param nDir Navigation direction prescription
  /// @param bcheck Boundary check for the navigation action
  /// @param sobject Start object to check against
  /// @param eobject End object to check against
  /// @param maxStepLength Maximal step length to check against
  NavigationOptions(NavigationDirection ndir, BoundaryCheck bcheck,
                    bool resolves = true, bool resolvem = true,
                    bool resolvep = false, const object_t* sobject = nullptr,
                    const object_t* eobject = nullptr)
      : navDir(ndir),
        boundaryCheck(std::move(bcheck)),
        resolveSensitive(resolves),
        resolveMaterial(resolvem),
        resolvePassive(resolvep),
        startObject(sobject),
        endObject(eobject),
        pathLimit(ndir * std::numeric_limits<double>::max()),
        overstepLimit(-1 * UnitConstants::um) {}
};
}  // namespace Acts
