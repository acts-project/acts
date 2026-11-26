// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts::detail {

/// Estimate the propagation direction to the closest straight-line surface
/// intersection, ignoring bounds. Zero path length selects forward propagation.
/// @param gctx Geometry context
/// @param surface Surface to approach
/// @param position Current global position
/// @param direction Current unit momentum direction
/// @return Forward or backward according to the signed intersection path length
inline Direction directionToSurface(const GeometryContext& gctx,
                                    const Surface& surface,
                                    const Vector3& position,
                                    const Vector3& direction) {
  const Intersection3D intersection =
      surface
          .intersect(gctx, position, direction, BoundaryTolerance::Infinite())
          .closest();
  return Direction::fromScalarZeroAsPositive(intersection.pathLength());
}

}  // namespace Acts::detail
