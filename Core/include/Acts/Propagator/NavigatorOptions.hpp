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
class TrackingVolume;

/// A surface of the tracking geometry with the bounds check relaxed to the
/// tolerance of the caller. It is targeted only inside its own volume (Gen3)
/// or on its layer (Gen1). To reach it outside, use an @c AdditionalSurface.
struct ExtendedSurface {
  /// The surface. The geometry must hold it, and it must outlive the
  /// propagation.
  const Surface* surface{};

  /// Tolerance used to intersect the surface. @c Infinite() drops the bounds
  /// check.
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::Infinite();
};

/// A surface the navigator offers on top of the tracking geometry, whether the
/// geometry holds it or not. It is intersected on every step and targeted
/// whenever it is closer than the next candidate of the geometry, which is held
/// back, not skipped. A closer portal wins, and a tie goes to the geometry.
struct AdditionalSurface {
  /// The surface. It must outlive the propagation.
  const Surface* surface{};

  /// Tolerance used to intersect the surface. @c Infinite() drops the bounds
  /// check.
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::Infinite();

  /// Volume the navigator offers the surface in. Null offers it everywhere.
  const TrackingVolume* volume{};

  /// Stop offering the surface once reached. Matters for a surface a track can
  /// approach twice, such as a line.
  bool dropAfterReached = true;
};

/// Plain navigator options carrying geometry context and navigation settings.
///
/// These are bound to the lifetime of a navigation state and have to be
/// invariant across all runs it serves. Per-run inputs belong into
/// @c NavigatorInitializeArguments instead.
struct NavigatorPlainOptions {
  /// NavigatorPlainOptions with context
  /// @param gctx The geometry context
  explicit NavigatorPlainOptions(const GeometryContext& gctx)
      : geoContext(gctx) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The surface tolerance
  double surfaceTolerance = s_onSurfaceTolerance;

  /// The near limit to resolve surfaces
  double nearLimit = s_onSurfaceTolerance;

  /// The far limit to resolve surfaces
  double farLimit = std::numeric_limits<double>::max();

  /// Surfaces of the tracking geometry with bounds extended by the caller
  std::vector<ExtendedSurface> extendedSurfaces;

  /// Relax the bounds check of a surface of the tracking geometry. By default
  /// the bounds check is dropped.
  /// @param surface The surface of the tracking geometry
  /// @param boundaryTolerance The tolerance used to intersect the surface
  void registerExtendedSurface(const Surface& surface,
                               const BoundaryTolerance& boundaryTolerance =
                                   BoundaryTolerance::Infinite()) {
    extendedSurfaces.emplace_back(&surface, boundaryTolerance);
  }

  /// Surfaces the navigator offers on top of the tracking geometry
  std::vector<AdditionalSurface> additionalSurfaces;

  /// Offer a surface on top of the tracking geometry. By default without a
  /// bounds check, in every volume, and only until reached.
  /// @param surface The surface to offer
  /// @param boundaryTolerance The tolerance used to intersect the surface
  /// @param volume The volume to offer the surface in, null for every volume
  /// @param dropAfterReached Whether to stop offering it once reached
  void registerAdditionalSurface(const Surface& surface,
                                 const BoundaryTolerance& boundaryTolerance =
                                     BoundaryTolerance::Infinite(),
                                 const TrackingVolume* volume = nullptr,
                                 bool dropAfterReached = true) {
    additionalSurfaces.emplace_back(&surface, boundaryTolerance, volume,
                                    dropAfterReached);
  }
};

}  // namespace Acts
