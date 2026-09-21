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

/// A surface of the tracking geometry whose bounds the navigator extends to the
/// tolerance of the caller. Only the bounds check changes, the geometry still
/// resolves the surface where it did before.
///
/// The extension applies only in the volume of the surface (Gen3), or on its
/// layer (Gen1). It can put the intersection outside that volume. The
/// navigator never targets such an intersection, because the propagation
/// leaves the volume through a boundary first, and no other volume offers the
/// surface. To reach a surface outside its volume, use an
/// @c AdditionalSurface instead.
struct ExtendedSurface {
  /// The surface. The geometry must hold it, and it must outlive the
  /// propagation.
  const Surface* surface{};

  /// Tolerance used to intersect the surface. @c Infinite() drops the bounds
  /// check.
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::Infinite();
};

/// A surface the navigator offers on top of the candidates of the tracking
/// geometry. The geometry can hold the surface or not.
///
/// The navigator targets it whenever it is closer than the candidate of the
/// geometry, and holds that candidate back rather than consuming it, so
/// nothing of the geometry is skipped. It intersects the surface on every
/// step, because @c Surface::intersect estimates along a straight line while
/// a track in a magnetic field bends away from it.
///
/// A boundary of the current volume is a candidate of the geometry too, so a
/// closer boundary wins. The propagation therefore reaches the surface only in
/// the volume that holds the intersection, and the current volume stays
/// correct.
///
/// If the geometry also offers the surface, a tie goes to the candidate of the
/// geometry, so the geometry keeps the handling of the surface in its own
/// volume.
struct AdditionalSurface {
  /// The surface. It must outlive the propagation.
  const Surface* surface{};

  /// Tolerance used to intersect the surface. @c Infinite() drops the bounds
  /// check.
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::Infinite();

  /// Volume the navigator offers the surface in. Null offers it everywhere.
  const TrackingVolume* volume{};

  /// Whether the navigator stops offering the surface once the propagation
  /// reached it. Only matters for a surface a track can approach twice, such
  /// as the point of closest approach to a line.
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

  /// Extend the bounds of a surface of the tracking geometry to the given
  /// tolerance. By default the bounds check is dropped. The navigator targets
  /// the surface only inside its volume, see @c ExtendedSurface.
  /// @param surface The surface of the tracking geometry
  /// @param boundaryTolerance The tolerance used to intersect the surface
  void registerExtendedSurface(const Surface& surface,
                               const BoundaryTolerance& boundaryTolerance =
                                   BoundaryTolerance::Infinite()) {
    extendedSurfaces.emplace_back(&surface, boundaryTolerance);
  }

  /// Surfaces the navigator offers on top of the tracking geometry
  std::vector<AdditionalSurface> additionalSurfaces;

  /// Offer a surface on top of the tracking geometry. By default the
  /// bounds check is dropped, every volume offers it, and the navigator stops
  /// offering it once the propagation reached it.
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
