// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

class RegularSurface : public Surface {
 public:
  using Surface::Surface;

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position where the normal vector is
  /// constructed
  ///
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector2& lposition) const = 0;

  Vector3 normal(const GeometryContext& gctx, const Vector2& pos,
                 const Vector3& /*direction*/) const override {
    return normal(gctx, pos);
  };

  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& /*direction*/,
      double tolerance = s_onSurfaceTolerance) const override {
    return globalToLocal(gctx, position, tolerance);
  }

  virtual Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      double tolerance = s_onSurfaceTolerance) const = 0;

  Vector3 localToGlobal(const GeometryContext& gctx, const Vector2& lposition,
                        const Vector3& /*direction*/) const override {
    return localToGlobal(gctx, lposition);
  }

  virtual Vector3 localToGlobal(const GeometryContext& gctx,
                                const Vector2& lposition) const = 0;
};

}  // namespace Acts
