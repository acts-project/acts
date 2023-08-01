// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RegularSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {

Vector3 RegularSurface::normal(const GeometryContext& gctx,
                               const Vector3& position) const {
  if (!isOnSurface(gctx, position, Vector3::Zero(), false)) {
    throw std::runtime_error("Not on surface");
  }
  return normal(gctx, Vector2{Vector2::Zero()});
}

Vector3 RegularSurface::normal(const GeometryContext& gctx) const {
  return normal(gctx, center(gctx));
}

Vector3 RegularSurface::normal(const GeometryContext& gctx, const Vector3& pos,
                               const Vector3& direction) const {
  throw_assert(isOnSurface(gctx, pos, direction, false), "Not on surface");
  return normal(gctx, pos);
};
}  // namespace Acts
