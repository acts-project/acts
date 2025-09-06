// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/BoundaryTolerance.hpp"

namespace Acts {

class Surface;

/// @brief The navigation target
///
/// This struct represents a navigation target which is communicated from the
/// navigator to the stepper through the propagator.
///
/// @note This incorporates `std::optional` semantics as the next target might
///       not exist.
struct NavigationTarget {
  /// Target surface for navigation
  const Surface* surface = nullptr;
  /// Index of the surface intersection to use
  std::uint8_t surfaceIntersectionIndex = 0;
  /// Boundary tolerance for intersection calculation
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

  /// Create an empty NavigationTarget (no target)
  /// @return Empty NavigationTarget indicating no navigation target
  static NavigationTarget None() { return NavigationTarget(); }

  /// Construct a NavigationTarget with a specific surface target.
  ///
  /// @param surface_ Target surface for navigation
  /// @param surfaceIntersectionIndex_ Index of the surface intersection to use
  /// @param boundaryTolerance_ Boundary tolerance for intersection calculation
  NavigationTarget(const Surface& surface_,
                   std::uint8_t surfaceIntersectionIndex_,
                   BoundaryTolerance boundaryTolerance_)
      : surface(&surface_),
        surfaceIntersectionIndex(surfaceIntersectionIndex_),
        boundaryTolerance(boundaryTolerance_) {}

  /// Check if this NavigationTarget represents no target
  /// @return True if no target is set, false otherwise
  bool isNone() const { return surface == nullptr; }

 private:
  NavigationTarget() = default;
};

}  // namespace Acts
