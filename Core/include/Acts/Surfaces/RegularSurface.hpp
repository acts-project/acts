// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

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

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param position is the global position where the normal vector is
  /// constructed
  /// @param gctx The current geometry context object, e.g. alignment

  ///
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector3& position) const;

  /// Return method for the normal vector of the surface
  ///
  /// It will return a normal vector at the center() position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  //
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx) const;

  Vector3 normal(const GeometryContext& gctx, const Vector3& pos,
                 const Vector3& direction) const override;
};
}  // namespace Acts
