// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <functional>

namespace Acts {

class GeometryContext;
class Surface;
class TrackingVolume;

/// Plain navigator options carrying geometry context and surfaces.
struct NavigatorPlainOptions {
  /// NavigatorPlainOptions with context
  /// @param gctx The geometry context
  explicit NavigatorPlainOptions(const GeometryContext& gctx)
      : geoContext(gctx) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// Start surface for navigation
  const Surface* startSurface{};
  /// Target surface for navigation
  const Surface* targetSurface{};

  /// The surface tolerance
  double surfaceTolerance = s_onSurfaceTolerance;

  /// The near limit to resolve surfaces
  double nearLimit = s_onSurfaceTolerance;

  /// The far limit to resolve surfaces
  double farLimit = std::numeric_limits<double>::max();

  /// Switch deciding whether surfaces without boundary
  /// are cleared after a volume switch although they've not been reached yet
  bool eraseUnboundVolChange{false};

  /// Delegate to decide whether free surfaces are appended to the navigation
  /// stream given the current volume and the track coordinates. If the
  /// delegate is set, it is called in each candidate resolution step
  /// for each surface that has not been marked as reached yet.
  /// @param gctx Current geometry context carrying the alignment information
  /// @param currentVol The current tracking volume in which the propagator resides
  /// @param pos Position of the track in global coordinates
  /// @param dir Direction vector of the track
  /// @param surface Free surface candidate to test
  using FreeSurfaceSelctor = Delegate<bool(
      const GeometryContext& gctx, const TrackingVolume& currentVol,
      const Vector3& pos, const Vector3& dir, const Surface& candidate)>;

  /// Delegate for selecting free surfaces during navigation
  FreeSurfaceSelctor freeSurfaceSelector;

  /// Surfaces that are not part of the tracking geometry
  std::vector<const Surface*> externalSurfaces;

  /// Append an external surface to be considered during navigation. The surface
  /// will be intersected without bounds check to force the propagation to
  /// resolve it. Note the surface must outlive the propagation and surfaces
  /// have to be added in the order they should be considered during navigation.
  /// @param surface The surface to add to the list
  void appendExternalSurface(const Surface& surface) {
    externalSurfaces.push_back(&surface);
  }
};

}  // namespace Acts
