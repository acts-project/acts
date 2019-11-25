// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlanarHelper.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

/// @brief Helpers for planar surfaces that share the same maths
namespace PlanarHelper {

/// Intersection with a planar surface
///
/// @param transform The 3D affine transform that places the surface
/// @param position The starting position for the intersection
/// @param direction The starting direction for the intersection
static Intersection intersectionEstimate(const Transform3D& transform,
                                         const Vector3D& position,
                                         const Vector3D& direction) {
  // Get the matrix from the transform (faster access)
  const auto& tMatrix = transform.matrix();
  const Vector3D pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3D pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // It is solvable, so go on
  double denom = direction.dot(pnormal);
  if (denom != 0.0) {
    // Translate that into a path
    double path = (pnormal.dot((pcenter - position))) / (denom);
    // Is valid hence either on surface or reachable
    Intersection::Status status =
        (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance)
            ? Intersection::Status::onSurface
            : Intersection::Status::reachable;
    // Return the intersection
    return Intersection{(position + path * direction), path, status};
  }
  return Intersection();
}

}  // namespace PlanarHelper
}  // namespace Acts
