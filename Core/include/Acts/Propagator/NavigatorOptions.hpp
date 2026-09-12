// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <functional>
#include <vector>

namespace Acts {

class GeometryContext;
class Surface;

/// A surface of the tracking geometry that the navigator intersects with the
/// tolerance of the caller. Only the bounds check changes, the geometry still
/// resolves the surface where it did before.
struct BoundaryToleranceOverride {
  /// The surface. The geometry must hold it, and it must outlive the
  /// propagation.
  const Surface* surface{};

  /// Tolerance used to intersect the surface. @c Infinite() drops the bounds
  /// check.
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::Infinite();
};

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

  /// Surfaces the navigator intersects with the tolerance of the caller
  std::vector<BoundaryToleranceOverride> boundaryToleranceOverrides;

  /// Intersect a surface of the tracking geometry with the given tolerance.
  /// By default the bounds check is dropped.
  /// @param surface The surface of the tracking geometry
  /// @param boundaryTolerance The tolerance used to intersect the surface
  void overrideBoundaryTolerance(const Surface& surface,
                                 const BoundaryTolerance& boundaryTolerance =
                                     BoundaryTolerance::Infinite()) {
    boundaryToleranceOverrides.push_back({&surface, boundaryTolerance});
  }
};

}  // namespace Acts
